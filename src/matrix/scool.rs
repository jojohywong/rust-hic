//! Build a single-cell single-resolution cooler (`.scool`) from CB pairs.
//!
//! Schema: <https://cooler.readthedocs.io/en/stable/schema.html#single-cell-single-resolution>
//!
//! `/chroms` and `/bins` are written once; every `/cells/<id>` group hard-links
//! them and stores its own `pixels` + `indexes`. Pixels are raw (unbalanced) and
//! only the upper triangle is stored.

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use hdf5::{
    types::{EnumMember, EnumType, FixedAscii, IntSize, TypeDescriptor},
    File as H5File, Group,
};
use rayon::prelude::*;
use std::{
    collections::{HashMap, HashSet},
    fs::File as StdFile,
    io::{BufRead, BufReader},
    path::Path,
    time::Instant,
};

use super::{
    attr::{int_attr, str_attr},
    read_chromosomes,
};

fn u64_attr(group: &Group, key: &str, value: u64) -> Result<()> {
    int_attr(group, key, i64::try_from(value).context("attribute exceeds i64")?)
}

fn group_name(cell: &str) -> String {
    let mut name = cell.replace(['/', '\0'], "_");
    if name.is_empty() || name == "." || name == ".." {
        name = "_cell".into();
    }
    name
}

fn write_i64(group: &Group, name: &str, data: &[i64]) -> Result<()> {
    if data.is_empty() {
        group.new_dataset_builder().with_data(data).create(name)?;
    } else {
        let chunk = data.len().min(65_536);
        group
            .new_dataset_builder()
            .with_data(data)
            .chunk((chunk,))
            .shuffle()
            .deflate(6)
            .create(name)?;
    }
    Ok(())
}

fn write_i32(group: &Group, name: &str, data: &[i32]) -> Result<()> {
    if data.is_empty() {
        group.new_dataset_builder().with_data(data).create(name)?;
    } else {
        let chunk = data.len().min(65_536);
        group
            .new_dataset_builder()
            .with_data(data)
            .chunk((chunk,))
            .shuffle()
            .deflate(6)
            .create(name)?;
    }
    Ok(())
}

fn open_pairs(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = StdFile::open(path).with_context(|| format!("open pairs {}", path.display()))?;
    if path.to_string_lossy().ends_with(".gz") {
        Ok(Box::new(BufReader::with_capacity(
            4 * 1024 * 1024,
            MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(4 * 1024 * 1024, file)))
    }
}

fn run() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 8 {
        bail!(
            "usage: build_scool PAIRS CHROMSIZES PER_CELL.tsv OUTPUT.scool BIN_SIZE MIN_PAIRS ASSEMBLY"
        )
    }
    let pairs = Path::new(&args[1]);
    let chromsizes = Path::new(&args[2]);
    let per_cell = Path::new(&args[3]);
    let output = Path::new(&args[4]);
    let bin_size: u64 = args[5].parse().context("invalid bin size")?;
    let min_pairs: u64 = args[6].parse().context("invalid min pairs")?;
    let assembly = &args[7];
    if bin_size == 0 {
        bail!("bin size must be positive")
    }
    let started = Instant::now();

    // Bins follow the same natural chromosome order as the mcool path.
    let chroms = read_chromosomes(chromsizes, true)?;
    let mut offsets: HashMap<&str, u64> = HashMap::with_capacity(chroms.len());
    let mut chrom_offset: Vec<i64> = Vec::with_capacity(chroms.len() + 1);
    let mut bin_chrom: Vec<i32> = Vec::new();
    let mut bin_start: Vec<i32> = Vec::new();
    let mut bin_end: Vec<i32> = Vec::new();
    let mut bins = 0u64;
    for (index, chrom) in chroms.iter().enumerate() {
        chrom_offset.push(i64::try_from(bins)?);
        offsets.insert(chrom.name.as_str(), bins);
        let mut start = 0u64;
        while start < chrom.len {
            let end = (start + bin_size).min(chrom.len);
            bin_chrom.push(i32::try_from(index)?);
            bin_start.push(i32::try_from(start)?);
            bin_end.push(i32::try_from(end)?);
            start = end;
            bins += 1;
        }
    }
    chrom_offset.push(i64::try_from(bins)?);
    if bins == 0 {
        bail!("no bins")
    }
    if bins > u32::MAX as u64 {
        bail!("too many bins for scool: {bins}")
    }
    let nchroms = chroms.len();

    // Cell set comes from the per-cell QC table; the order is deterministic.
    let mut cells: Vec<String> = Vec::new();
    let mut cell_index: HashMap<String, usize> = HashMap::new();
    for line in BufReader::new(StdFile::open(per_cell)?).lines() {
        let line = line?;
        if line.starts_with("cell_id") || line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let id = fields.next().unwrap_or("");
        let pairs: u64 = fields.next().unwrap_or("0").parse().unwrap_or(0);
        if id == "unassigned" || id.is_empty() || pairs < min_pairs {
            continue;
        }
        if !cell_index.contains_key(id) {
            cell_index.insert(id.to_owned(), cells.len());
            cells.push(id.to_owned());
        }
    }
    if cells.is_empty() {
        bail!("no cells pass --min-pairs {min_pairs}")
    }

    // Read each pair once. The input is position-sorted, so each cell's packed
    // pixel keys only need a final in-memory sort before run-length merging.
    let mut pixels: Vec<Vec<u64>> = vec![Vec::new(); cells.len()];
    let mut selected = 0u64;
    for line in open_pairs(pairs)?.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let mut f = [""; 9];
        let mut n = 0;
        for field in line.split('\t') {
            if n == 9 {
                break;
            }
            f[n] = field;
            n += 1;
        }
        if n != 9 {
            continue;
        }
        let Some(&cell) = cell_index.get(f[8]) else {
            continue;
        };
        let (Some(&off1), Some(&off2)) = (offsets.get(f[1]), offsets.get(f[3])) else {
            continue;
        };
        let (Ok(p1), Ok(p2)) = (f[2].parse::<u64>(), f[4].parse::<u64>()) else {
            continue;
        };
        if p1 == 0 || p2 == 0 {
            continue;
        }
        let mut b1 = off1 + (p1 - 1) / bin_size;
        let mut b2 = off2 + (p2 - 1) / bin_size;
        if b1 > b2 {
            std::mem::swap(&mut b1, &mut b2);
        }
        if b1 >= bins || b2 >= bins {
            continue;
        }
        pixels[cell].push((b1 << 32) | b2);
        selected += 1;
    }
    eprintln!(
        "event=scool_scan cells={} pairs={selected} bins={bins} elapsed_seconds={:.3}",
        cells.len(),
        started.elapsed().as_secs_f64()
    );

    pixels.par_iter_mut().for_each(|v| v.par_sort_unstable());

    let file = H5File::create(output)
        .with_context(|| format!("create scool {}", output.display()))?;
    let root = file.group("/")?;
    let chrom_group = root.create_group("chroms")?;
    let bin_group = root.create_group("bins")?;
    let cells_group = root.create_group("cells")?;

    let names: Vec<FixedAscii<256>> = chroms
        .iter()
        .map(|c| FixedAscii::<256>::from_ascii(&c.name))
        .collect::<std::result::Result<_, _>>()?;
    let lengths: Vec<i32> = chroms
        .iter()
        .map(|c| i32::try_from(c.len))
        .collect::<std::result::Result<_, _>>()?;
    let chrom_chunk = names.len().max(1);
    chrom_group
        .new_dataset_builder()
        .with_data(&names)
        .chunk((chrom_chunk,))
        .shuffle()
        .deflate(6)
        .create("name")?;
    chrom_group
        .new_dataset_builder()
        .with_data(&lengths)
        .chunk((chrom_chunk,))
        .shuffle()
        .deflate(6)
        .create("length")?;

    let chrom_type = TypeDescriptor::Enum(EnumType {
        size: IntSize::U4,
        signed: true,
        members: chroms
            .iter()
            .enumerate()
            .map(|(i, c)| EnumMember {
                name: c.name.clone(),
                value: i as u64,
            })
            .collect(),
    });
    let chrom_ds = bin_group
        .new_dataset_builder()
        .empty_as(&chrom_type)
        .shape(bin_chrom.len())
        .chunk((65_536.min(bin_chrom.len().max(1)),))
        .shuffle()
        .deflate(6)
        .create("chrom")?;
    let chrom_dtype = chrom_ds.dtype()?;
    let status = hdf5::sync::sync(|| unsafe {
        hdf5_sys::h5d::H5Dwrite(
            chrom_ds.id(),
            chrom_dtype.id(),
            hdf5_sys::h5s::H5S_ALL,
            hdf5_sys::h5s::H5S_ALL,
            hdf5_sys::h5p::H5P_DEFAULT,
            bin_chrom.as_ptr().cast(),
        )
    });
    if status < 0 {
        bail!("failed to write chromosome enum values")
    }
    let bin_chunk = bin_start.len().max(1).min(65_536);
    bin_group
        .new_dataset_builder()
        .with_data(&bin_start)
        .chunk((bin_chunk,))
        .shuffle()
        .deflate(6)
        .create("start")?;
    bin_group
        .new_dataset_builder()
        .with_data(&bin_end)
        .chunk((bin_chunk,))
        .shuffle()
        .deflate(6)
        .create("end")?;
    let weight = vec![1.0f64; bin_start.len()];
    bin_group
        .new_dataset_builder()
        .with_data(&weight)
        .chunk((bin_chunk,))
        .shuffle()
        .deflate(6)
        .create("weight")?;

    str_attr(&root, "format", "HDF5::SCOOL")?;
    int_attr(&root, "format-version", 1)?;
    str_attr(&root, "bin-type", "fixed")?;
    u64_attr(&root, "bin-size", bin_size)?;
    str_attr(&root, "storage-mode", "symmetric-upper")?;
    str_attr(&root, "assembly", assembly)?;
    str_attr(&root, "generated-by", "rust-hic")?;
    u64_attr(&root, "nbins", bins)?;
    u64_attr(&root, "nchroms", nchroms as u64)?;
    u64_attr(&root, "ncells", cells.len() as u64)?;

    let mut used_names = HashSet::new();
    let mut total_pixels = 0u64;
    for (cell, name) in cells.iter().enumerate() {
        let mut gname = group_name(name);
        if !used_names.insert(gname.clone()) {
            let mut k = 2usize;
            loop {
                let candidate = format!("{gname}__{k}");
                if used_names.insert(candidate.clone()) {
                    gname = candidate;
                    break;
                }
                k += 1;
            }
        }
        let group = cells_group.create_group(&gname)?;
        group.link_hard("/bins", "bins")?;
        group.link_hard("/chroms", "chroms")?;

        let keys = &pixels[cell];
        let mut b1v: Vec<i64> = Vec::new();
        let mut b2v: Vec<i64> = Vec::new();
        let mut counts: Vec<i32> = Vec::new();
        let mut bin_offset: Vec<i64> = Vec::with_capacity(bins as usize + 1);
        bin_offset.push(0);
        let mut cursor = 0u64;
        let mut nnz = 0i64;
        let mut i = 0usize;
        while i < keys.len() {
            let key = keys[i];
            let mut j = i + 1;
            while j < keys.len() && keys[j] == key {
                j += 1;
            }
            let b1 = key >> 32;
            while cursor < b1 {
                bin_offset.push(nnz);
                cursor += 1;
            }
            b1v.push(i64::try_from(b1)?);
            b2v.push(i64::try_from(key & 0xffff_ffff)?);
            counts.push(i32::try_from(j - i).context("pixel count exceeds i32")?);
            nnz += 1;
            i = j;
        }
        while cursor < bins {
            bin_offset.push(nnz);
            cursor += 1;
        }
        total_pixels += nnz as u64;

        let pixel_group = group.create_group("pixels")?;
        write_i64(&pixel_group, "bin1_id", &b1v)?;
        write_i64(&pixel_group, "bin2_id", &b2v)?;
        write_i32(&pixel_group, "count", &counts)?;
        let index_group = group.create_group("indexes")?;
        write_i64(&index_group, "bin1_offset", &bin_offset)?;
        write_i64(&index_group, "chrom_offset", &chrom_offset)?;

        str_attr(&group, "format", "HDF5::Cooler")?;
        int_attr(&group, "format-version", 3)?;
        str_attr(&group, "bin-type", "fixed")?;
        u64_attr(&group, "bin-size", bin_size)?;
        str_attr(&group, "storage-mode", "symmetric-upper")?;
        str_attr(&group, "assembly", assembly)?;
        u64_attr(&group, "nbins", bins)?;
        u64_attr(&group, "nchroms", nchroms as u64)?;
        u64_attr(&group, "nnz", nnz as u64)?;
    }

    eprintln!(
        "event=scool_done cells={} pairs={selected} pixels={total_pixels} bins={bins} path={} elapsed_seconds={:.3}",
        cells.len(),
        output.display(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

pub(crate) fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error:#}");
        std::process::exit(1)
    }
}
