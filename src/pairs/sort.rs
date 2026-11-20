use super::format::Format;
use super::Key;
use anyhow::{bail, Context, Result};
use std::{
    collections::HashMap,
    env,
    fs::{self, File, OpenOptions},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::Stdio,
    thread,
    time::Instant,
};
struct Row {
    key: Key,
    line: String,
}
fn write_run(mut rows: Vec<Row>, path: PathBuf, header: String) -> Result<()> {
    rows.sort_by(|a, b| a.key.cmp(&b.key));
    let f = File::create(path)?;
    let mut out = BufWriter::with_capacity(1024 * 1024, f);
    out.write_all(header.as_bytes())?;
    for row in rows {
        out.write_all(row.line.as_bytes())?
    }
    out.flush()?;
    Ok(())
}
pub(crate) fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
fn canonical_rank(name: &str) -> (u8, u32, u32) {
    let short = name.strip_prefix("chr").unwrap_or(name);
    match short {
        "X" => (1, 0, 0),
        "Y" => (2, 0, 0),
        "Z" => (3, 0, 0),
        "W" => (4, 0, 0),
        "M" | "MT" => (5, 0, 0),
        value => {
            let digits = value.bytes().take_while(|b| b.is_ascii_digit()).count();
            if digits == 0 {
                return (6, 0, 0);
            }
            let number = value[..digits].parse::<u32>().unwrap_or(0);
            if number == 0 {
                return (6, 0, 0);
            }
            // Drosophila-style arms: 2L, 2R, 3L, 3R come before 4/X/Y.
            let sub = match &value[digits..] {
                "" => 0,
                "L" => 1,
                "R" => 2,
                suffix => 3 + u32::from(suffix.as_bytes()[0]),
            };
            (0, number, sub)
        }
    }
}

fn natural_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    let mut left = a.as_bytes();
    let mut right = b.as_bytes();
    while !left.is_empty() && !right.is_empty() {
        if left[0].is_ascii_digit() && right[0].is_ascii_digit() {
            let left_end = left.iter().position(|c| !c.is_ascii_digit()).unwrap_or(left.len());
            let right_end = right.iter().position(|c| !c.is_ascii_digit()).unwrap_or(right.len());
            let left_number: u128 = std::str::from_utf8(&left[..left_end]).unwrap_or("0").parse().unwrap_or(0);
            let right_number: u128 = std::str::from_utf8(&right[..right_end]).unwrap_or("0").parse().unwrap_or(0);
            match left_number.cmp(&right_number) {
                std::cmp::Ordering::Equal => {}
                other => return other,
            }
            left = &left[left_end..];
            right = &right[right_end..];
        } else {
            match left[0].cmp(&right[0]) {
                std::cmp::Ordering::Equal => {}
                other => return other,
            }
            left = &left[1..];
            right = &right[1..];
        }
    }
    left.len().cmp(&right.len())
}

#[cfg(test)]
mod tests {
    use super::{canonical_rank, natural_cmp};

    #[test]
    fn chromosome_order_is_natural() {
        assert!(canonical_rank("chr2") < canonical_rank("chr10"));
        assert!(canonical_rank("chr10") < canonical_rank("chrX"));
        assert!(canonical_rank("chrX") < canonical_rank("chrM"));
        // Drosophila melanogaster arms.
        assert!(canonical_rank("chr2L") < canonical_rank("chr2R"));
        assert!(canonical_rank("chr2R") < canonical_rank("chr3L"));
        assert!(canonical_rank("chr3R") < canonical_rank("chr4"));
        assert!(canonical_rank("chr4") < canonical_rank("chrX"));
        assert_eq!(
            natural_cmp("chr2_random", "chr10_random"),
            std::cmp::Ordering::Less
        );
    }
}

fn run() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 8 {
        bail!("usage: sort_pairs SIZES ASSEMBLY OUTPUT TMPDIR THREADS RUN_MIB INPUT.pairs")
    };
    let start = Instant::now();
    let output = Path::new(&args[3]);
    if output.exists() {
        bail!("output exists")
    };
    let threads: usize = args[5].parse()?;
    let limit: usize = args[6]
        .parse::<usize>()?
        .checked_mul(1024 * 1024)
        .context("run size overflow")?;
    if threads == 0 || limit < 1024 * 1024 {
        bail!("invalid thread count or run MiB")
    };
    let mut chroms: Vec<(String, u64)> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for line in BufReader::new(File::open(&args[1])?).lines() {
        let l = line?;
        let a: Vec<_> = l.split_whitespace().collect();
        if a.is_empty() {
            continue;
        }
        if a.len() != 2 {
            bail!("invalid chromsizes")
        };
        let n: u64 = a[1].parse()?;
        if n == 0 || !seen.insert(a[0].to_owned()) {
            bail!("invalid chromsize")
        };
        chroms.push((a[0].to_owned(), n))
    }
    if chroms.is_empty() {
        bail!("empty chromsizes")
    };
    // Natural chromosome order (chr1, chr2, ..., chr10, chrX, chrY, chrM) in
    // every mode. Pairtools comparisons must use a sizes file in this same
    // order for exact tie-breaking parity.
    let original_order: Vec<String> = chroms.iter().map(|(name, _)| name.clone()).collect();
    chroms.sort_by(|a, b| {
        canonical_rank(&a.0)
            .cmp(&canonical_rank(&b.0))
            .then_with(|| natural_cmp(&a.0, &b.0))
    });
    let natural_order: Vec<String> = chroms.iter().map(|(name, _)| name.clone()).collect();
    if natural_order != original_order {
        eprintln!(
            "event=chromsizes_reordered mode=natural original={} natural={}",
            original_order.join(","),
            natural_order.join(",")
        );
    }
    let mut order = HashMap::new();
    for (index, (name, _)) in chroms.iter().enumerate() {
        order.insert(name.clone(), index);
    }
    if args[2].contains(['\r', '\n']) {
        bail!("invalid assembly")
    };
    let mut header=format!("## pairs format v1.0.0\n#sorted: chr1-chr2-pos1-pos2\n#shape: upper triangle\n#genome_assembly: {}\n",args[2]);
    for (c, n) in &chroms {
        header += &format!("#chromsize: {c} {n}\n")
    }
    let scratch = Path::new(&args[4]).join(format!("rust-sort-{}", std::process::id()));
    fs::create_dir(&scratch).context("create exclusive scratch directory")?;
    let mut temp = output.to_path_buf();
    let name = output
        .file_name()
        .context("missing output filename")?
        .to_string_lossy();
    temp.set_file_name(format!(".{name}.{}.partial", std::process::id()));
    let result = (|| -> Result<()> {
        let mut jobs: Vec<thread::JoinHandle<Result<()>>> = Vec::new();
        let mut runs = Vec::new();
        let mut rows = Vec::new();
        let mut size = 0usize;
        let mut total = 0u64;
        let (mut body, mut columns) = (false, None);
        let mut reader = BufReader::with_capacity(1024 * 1024, File::open(&args[7])?);
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if line.starts_with('#') {
                if body {
                    bail!("header after body")
                };
                if let Some(meta) = line.strip_prefix("#chromsize:") {
                    let parts: Vec<_> = meta.split_whitespace().collect();
                    if parts.len() != 2 {
                        bail!("invalid #chromsize header");
                    }
                    let observed: u64 = parts[1].parse().context("invalid #chromsize length")?;
                    let idx = *order.get(parts[0]).with_context(|| {
                        format!("header chromosome absent from sizes: {}", parts[0])
                    })?;
                    if chroms[idx].1 != observed {
                        bail!("header length disagrees with reference for {}", parts[0]);
                    }
                }
                if line.starts_with("#columns:") {
                    if columns.is_some() {
                        bail!("duplicate #columns header")
                    }
                    let format = Format::from_columns(&line)?;
                    header.push_str(&format.columns_line());
                    columns = Some(format);
                }
                continue;
            }
            if line.trim().is_empty() {
                continue;
            }
            body = true;
            let format = columns.context("missing #columns header")?;
            let v = format.fields(&line)?;
            let p1: u64 = v[2].parse()?;
            let p2: u64 = v[4].parse()?;
            if (v[1] != "!" && !order.contains_key(v[1]))
                || (v[3] != "!" && !order.contains_key(v[3]))
            {
                bail!("unknown chromosome {} or {}", v[1], v[3])
            }
            let a = order.get(v[1]).copied();
            let b = order.get(v[3]).copied();
            let known = a.is_some() && b.is_some() && p1 > 0 && p2 > 0;
            if known {
                let (ia, ib) = (a.unwrap(), b.unwrap());
                if p1 > chroms[ia].1 || p2 > chroms[ib].1 {
                    bail!("invalid pair coordinates")
                }
            }
            if !["+", "-"].contains(&v[5]) || !["+", "-"].contains(&v[6]) {
                bail!("invalid pair strands")
            };
            // Unmapped (`!`) and contigs absent from the sizes file are kept and
            // sorted without flipping, matching pairtools.
            let flip = if known {
                (a.unwrap(), p1) > (b.unwrap(), p2)
            } else {
                false
            };
            let (c1, x, c2, y, s1, s2) = if flip {
                (v[3], p2, v[1], p1, v[6], v[5])
            } else {
                (v[1], p1, v[3], p2, v[5], v[6])
            };
            let kind = v[7];
            if kind.len() != 2 || !kind.is_ascii() {
                bail!("pair_type must contain two ASCII characters")
            }
            let kind = if flip {
                kind.chars().rev().collect()
            } else {
                kind.to_owned()
            };
            let mut normalized = format!("{}\t{c1}\t{x}\t{c2}\t{y}\t{s1}\t{s2}\t{kind}", v[0]);
            if format == Format::Fragments {
                let (left, right) = if flip {
                    (&v[11..14], &v[8..11])
                } else {
                    (&v[8..11], &v[11..14])
                };
                for field in left.iter().chain(right) {
                    normalized.push('\t');
                    normalized.push_str(field);
                }
            } else if format == Format::Cell {
                normalized.push('\t');
                normalized.push_str(v[8]);
            }
            normalized.push('\n');
            size += normalized.len();
            rows.push(Row {
                key: Key {
                    c1: c1.into(),
                    c2: c2.into(),
                    p1: x,
                    p2: y,
                    kind,
                },
                line: normalized,
            });
            total += 1;
            if size >= limit {
                let path = scratch.join(format!("run-{:06}.pairs", runs.len()));
                runs.push(path.clone());
                eprintln!(
                    "event=run_started run={} records={total} elapsed_seconds={:.3}",
                    runs.len(),
                    start.elapsed().as_secs_f64()
                );
                let batch = std::mem::take(&mut rows);
                let h = header.clone();
                jobs.push(thread::spawn(move || write_run(batch, path, h)));
                size = 0;
                if jobs.len() >= threads {
                    jobs.remove(0)
                        .join()
                        .map_err(|_| anyhow::anyhow!("sort worker panicked"))??
                }
            }
        }
        if columns.is_none() {
            bail!("missing columns header")
        };
        if !rows.is_empty() || runs.is_empty() {
            let path = scratch.join(format!("run-{:06}.pairs", runs.len()));
            runs.push(path.clone());
            let h = header.clone();
            jobs.push(thread::spawn(move || write_run(rows, path, h)))
        }
        for job in jobs {
            job.join()
                .map_err(|_| anyhow::anyhow!("sort worker panicked"))??
        }
        eprintln!(
            "event=runs_done records={total} runs={} elapsed_seconds={:.3}",
            runs.len(),
            start.elapsed().as_secs_f64()
        );
        let file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&temp)?;
        let status = crate::stage::stage_command(crate::stage::Stage::MergeRuns)?
            .args(&runs)
            .stdout(Stdio::from(file))
            .status()?;
        if !status.success() {
            bail!("run merge failed {status}")
        };
        fs::hard_link(&temp, output).context("publish output")?;
        eprintln!(
            "event=sort_done records={total} elapsed_seconds={:.3}",
            start.elapsed().as_secs_f64()
        );
        Ok(())
    })();
    if result.is_ok() {
        fs::remove_file(&temp)?;
        fs::remove_dir_all(&scratch).context("remove owned scratch directory")?;
    } else {
        eprintln!(
            "event=sort_failed scratch={} partial={}",
            scratch.display(),
            temp.display()
        );
    }
    result
}
