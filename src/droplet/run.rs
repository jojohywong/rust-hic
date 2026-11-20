//! Droplet Hi-C pipeline: index R1/R3, map R2 once, then align chunks.
use anyhow::{bail, Context, Result};
use std::{
    collections::HashMap,
    ffi::CString,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    os::unix::ffi::OsStrExt,
    os::unix::net::UnixStream,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Mutex,
    },
    thread,
    time::Instant,
};

use crate::{
    fastq_index::stream_split,
    mapping::{bam, sam},
    stage::{stage_command, Stage},
};

use super::{
    barcode,
    config::DropletConfig,
    fragments,
    index::{self, IndexManifest},
    qc,
};

const PAIRS_HEADER: &str =
    "## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type CB\n";

pub(crate) struct ChunkResult {
    pub(crate) pairs: PathBuf,
    pub(crate) bam: Option<PathBuf>,
    pub(crate) input_pairs: u64,
    pub(crate) selected_pairs: u64,
    pub(crate) barcode_missing: u64,
}

pub(crate) fn run(
    path: &Path,
    check: bool,
    keep_bam_flag: bool,
    barcode_preprocessing_only: bool,
) -> Result<()> {
    let (c, _base) = DropletConfig::load(path)?;
    if check {
        println!("Droplet configuration valid");
        println!("  R1: {}", c.read1.display());
        println!("  R2 (barcode): {}", c.barcode_read.display());
        println!("  R3: {}", c.read3.display());
        println!("  Reference: {}", c.reference.display());
        println!("  Index: {}", c.aligner_index.display());
        println!("  Aligner: {}", c.aligner_executable.display());
        println!("  fastp: {}", c.fastp_executable.display());
        println!("  Whitelist: {}", c.whitelist.display());
        println!(
            "  Barcode: slice {}:{} orientation {} max_mismatch {}",
            c.barcode.slice_start,
            c.barcode.slice_start + c.barcode.slice_len,
            c.barcode.orientation,
            c.barcode.max_mismatch
        );
        println!("  Output: {}", c.output_dir.display());
        println!("  Chunks: target {}", c.chunk_size);
        println!(
            "  Concurrency: {} jobs x {} align threads",
            c.jobs, c.threads_per_job
        );
        if c.index().join("manifest.json").is_file() {
            let m = IndexManifest::load(&c.index())?;
            println!("  Existing index: {} chunks", m.chunks);
        }
        return Ok(());
    }
    let started = Instant::now();
    let logs = c.output_dir.join("logs");
    let qc_dir = c.output_dir.join("qc");
    create_dirs(&c)?;

    let index_dir = index::ensure_index(&c, &logs.join("fastq-index.log"))?;
    let manifest = IndexManifest::load(&index_dir)?;
    let sidecar_dir = c.scratch().join("barcodes");
    if barcode_preprocessing_only {
        let preprocessing = barcode::preprocess(&c, &manifest, &sidecar_dir)?;
        write_barcode_qc(&preprocessing, &qc_dir)?;
        return Ok(());
    }
    let keep_bam = c.keep_bam || keep_bam_flag;
    let chunk_dir = c.scratch().join("chunks");
    fs::create_dir_all(&chunk_dir)?;

    let post = c.cell_barcode_process_mode == "post";
    let none = c.cell_barcode_process_mode == "none";
    let results = if none {
        // No barcode read at all: align only, leave CB as ".".
        let align_start = Instant::now();
        let results = align_chunks(
            &c,
            &index_dir,
            &sidecar_dir,
            manifest.chunks,
            false,
            keep_bam,
            &chunk_dir,
            &logs,
        )?;
        eprintln!(
            "event=align_done elapsed_seconds={:.3}",
            align_start.elapsed().as_secs_f64()
        );
        results
    } else if post && c.barcode_preprocessing_overlap {
        // Run the preprocessing in a child process so rapidgzip and the aligner never
        // share an address space (the in-process version crashed with SIGSEGV).
        let mut child = Command::new(std::env::current_exe()?)
            .arg("droplet")
            .arg(path)
            .arg("--barcode-preprocessing-only")
            .spawn()
            .context("start barcode preprocessing child")?;
        let align_start = Instant::now();
        let results = align_chunks(
            &c,
            &index_dir,
            &sidecar_dir,
            manifest.chunks,
            false,
            keep_bam,
            &chunk_dir,
            &logs,
        )?;
        eprintln!(
            "event=align_done elapsed_seconds={:.3}",
            align_start.elapsed().as_secs_f64()
        );
        let status = child.wait().context("wait for barcode preprocessing child")?;
        if !status.success() {
            bail!("barcode preprocessing child failed: {status}")
        }
        results
    } else {
        let preprocessing = barcode::preprocess(&c, &manifest, &sidecar_dir)?;
        write_barcode_qc(&preprocessing, &qc_dir)?;
        let align_start = Instant::now();
        let results = align_chunks(
            &c,
            &index_dir,
            &sidecar_dir,
            manifest.chunks,
            !post,
            keep_bam,
            &chunk_dir,
            &logs,
        )?;
        eprintln!(
            "event=align_done elapsed_seconds={:.3}",
            align_start.elapsed().as_secs_f64()
        );
        results
    };
    finish_pipeline(&c, started, &results, keep_bam)
}

/// Create the run directories and record the resolved configuration.
pub(crate) fn create_dirs(c: &DropletConfig) -> Result<()> {
    for dir in [
        &c.output_dir,
        &c.output_dir.join("logs"),
        &c.output_dir.join("qc"),
        &c.output_dir.join("alignments"),
        &c.output_dir.join("pairs"),
    ] {
        fs::create_dir_all(dir).with_context(|| format!("create {}", dir.display()))?;
    }
    serde_yaml::to_writer(File::create(c.output_dir.join("config.resolved.yaml"))?, c)?;
    Ok(())
}

/// Join post-alignment CB columns, merge BAMs, deduplicate, QC and build the matrix.
pub(crate) fn finish_pipeline(
    c: &DropletConfig,
    started: Instant,
    results: &[ChunkResult],
    keep_bam: bool,
) -> Result<()> {
    let logs = c.output_dir.join("logs");
    let qc_dir = c.output_dir.join("qc");
    let align_dir = c.output_dir.join("alignments");
    let pairs_dir = c.output_dir.join("pairs");
    if c.cell_barcode_process_mode == "post" {
        let join_start = Instant::now();
        join_chunks(results, &c.scratch().join("barcodes"), c.jobs)?;
        eprintln!(
            "event=cb_join_done elapsed_seconds={:.3}",
            join_start.elapsed().as_secs_f64()
        );
    }

    if keep_bam {
        let parts: Vec<_> = results
            .iter()
            .filter_map(|r| r.bam.as_ref().cloned())
            .collect();
        let out = align_dir.join(format!("{}.bam", c.prefix));
        bam::merge(&parts, &out).with_context(|| format!("merge BAM to {}", out.display()))?;
        eprintln!(
            "event=bam_merged parts={} path={}",
            parts.len(),
            out.display()
        );
        if c.atac_fragments {
            let fragments = c
                .output_dir
                .join("atac")
                .join(format!("{}.fragments.tsv.gz", c.prefix));
            let n = fragments::write(
                &out,
                &fragments,
                &fragments::FragmentConfig {
                    t5_shift: c.t5_shift,
                    min_length: c.atac_min_length,
                    max_length: c.atac_max_length,
                },
            )?;
            eprintln!(
                "event=fragments_done fragments={} path={}",
                n,
                fragments.display()
            );
        }
    } else if c.atac_fragments {
        bail!("atac_fragments requires keep_bam")
    }

    let mut pairs_paths = Vec::with_capacity(results.len());
    let (mut input_pairs, mut selected_pairs, mut barcode_missing) = (0u64, 0u64, 0u64);
    for r in results {
        pairs_paths.push(r.pairs.clone());
        input_pairs += r.input_pairs;
        selected_pairs += r.selected_pairs;
        barcode_missing += r.barcode_missing;
    }
    let merged = c.scratch().join("merged.pairs");
    merge_pairs(&pairs_paths, &merged)?;
    let sorted = c.scratch().join("sorted.pairs");
    sort_pairs(c, &merged, &sorted, &logs)?;
    // `keep_pairs` says the final pairs are a deliverable; `compressed_pairs`
    // selects the compressed form. When compression is on it supersedes the
    // plain file, so `keep_pairs: true` + `compressed_pairs: true` retains only
    // `<prefix>.pairs.gz` (+ `.px2`).
    let retain_pairs = c.keep_pairs || c.compressed_pairs;
    let keep_plain = c.keep_pairs && !c.compressed_pairs;
    let final_pairs = if retain_pairs {
        pairs_dir.join(format!("{}.pairs", c.prefix))
    } else {
        c.scratch().join("final.pairs")
    };
    dedup_pairs(c, &sorted, &final_pairs, &logs)?;

    let summary = qc::per_cell(&final_pairs, &qc_dir)?;
    if c.matrix.enabled {
        build_matrix(c, &final_pairs, &logs)?;
    }
    let mut pairs_gz = None;
    if c.compressed_pairs {
        let gz = pairs_dir.join(format!("{}.pairs.gz", c.prefix));
        crate::workflow::pairix::pack(
            &final_pairs,
            &gz,
            c.dedup.workers.clamp(1, 64),
            Some(&logs.join("pairix.log")),
        )?;
        pairs_gz = Some(gz);
    }
    write_pipeline_qc(
        &qc_dir,
        &summary,
        input_pairs,
        selected_pairs,
        barcode_missing,
        started.elapsed().as_secs_f64(),
    )?;
    write_manifest(
        c,
        keep_plain.then_some(final_pairs.as_path()),
        pairs_gz.as_deref(),
        &summary,
    )?;
    fs::write(
        c.output_dir.join("status.json"),
        format!(
            "{{\"status\":\"complete\",\"elapsed_seconds\":{:.3},\"pairs\":{}}}\n",
            started.elapsed().as_secs_f64(),
            summary.total_pairs
        ),
    )?;
    if !keep_plain && final_pairs.exists() {
        fs::remove_file(&final_pairs)
            .with_context(|| format!("remove unretained pairs {}", final_pairs.display()))?;
    }
    eprintln!(
        "event=droplet_done pairs={} cells={} elapsed_seconds={:.1}",
        summary.total_pairs,
        summary.cells,
        started.elapsed().as_secs_f64()
    );
    Ok(())
}


pub(crate) fn write_barcode_qc(preprocessing: &barcode::BarcodePreprocessing, qc_dir: &Path) -> Result<()> {
    let mut out = BufWriter::new(File::create(qc_dir.join("barcode_stats.tsv"))?);
    writeln!(out, "#orientation: {}", preprocessing.orientation)?;
    writeln!(
        out,
        "chunk\ttotal\texact_match\tmatch_within_threshold\tambiguous\tunassigned\tassignable_fraction"
    )?;
    for (i, s) in preprocessing.stats.iter().enumerate() {
        let assignable = if s.total == 0 {
            0.0
        } else {
            (s.exact + s.corrected) as f64 / s.total as f64
        };
        writeln!(
            out,
            "{i}\t{}\t{}\t{}\t{}\t{}\t{assignable:.6}",
            s.total,
            s.exact,
            s.corrected,
            s.ambiguous,
            s.unassigned + s.poly_g
        )?;
    }
    let t = &preprocessing.total;
    let assignable = if t.total == 0 {
        0.0
    } else {
        (t.exact + t.corrected) as f64 / t.total as f64
    };
    writeln!(
        out,
        "total\t{}\t{}\t{}\t{}\t{}\t{assignable:.6}",
        t.total,
        t.exact,
        t.corrected,
        t.ambiguous,
        t.unassigned + t.poly_g
    )?;
    out.flush()?;
    Ok(())
}

fn fifo(path: &Path) -> Result<()> {
    let name = CString::new(path.as_os_str().as_bytes())?;
    if unsafe { libc::mkfifo(name.as_ptr(), 0o600) } != 0 {
        return Err(std::io::Error::last_os_error()).context("create FIFO");
    }
    Ok(())
}

fn load_barcodes(path: &Path) -> Result<HashMap<Vec<u8>, Vec<u8>>> {
    let mut map = HashMap::new();
    let mut reader = BufReader::with_capacity(1024 * 1024, File::open(path)?);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let end = line.len() - usize::from(line.last() == Some(&b'\n'));
        let Some(tab) = line[..end].iter().position(|b| *b == b'\t') else {
            bail!("invalid barcode sidecar line in {}", path.display())
        };
        map.insert(line[..tab].to_vec(), line[tab + 1..end].to_vec());
    }
    Ok(map)
}

fn join_chunks(results: &[ChunkResult], sidecar_dir: &Path, workers: usize) -> Result<()> {
    if results.is_empty() {
        return Ok(());
    }
    let next = AtomicUsize::new(0);
    let errors = Mutex::new(Vec::new());
    thread::scope(|scope| {
        for _ in 0..workers.min(results.len()).max(1) {
            let next = &next;
            let errors = &errors;
            scope.spawn(move || loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= results.len() {
                    break;
                }
                let sidecar = sidecar_dir.join(format!("chunk-{i:04}.tsv"));
                if let Err(e) = join_pairs_cb(&results[i].pairs, &sidecar) {
                    errors.lock().unwrap().push(format!("chunk {i}: {e:#}"));
                }
            });
        }
    });
    let errors = errors.into_inner().unwrap();
    if !errors.is_empty() {
        bail!("{}", errors.join("; "))
    }
    Ok(())
}

/// Replace the barcode columns of an unstamped pairs file from its sidecar.
fn join_pairs_cb(pairs: &Path, sidecar: &Path) -> Result<()> {
    let map = load_barcodes(sidecar)?;
    let tmp = pairs.with_extension("cbjoin");
    let mut reader = BufReader::with_capacity(1024 * 1024, File::open(pairs)?);
    let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(&tmp)?);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.starts_with(b"#") {
            writer.write_all(&line)?;
            continue;
        }
        let mut fields = line.splitn(9, |b| *b == b'\t');
        let read = fields.next().unwrap_or(b"");
        let mut rest: [&[u8]; 7] = [b""; 7];
        for slot in rest.iter_mut() {
            *slot = fields.next().unwrap_or(b"");
        }
        let cb = map.get(read).map(|v| v.as_slice()).unwrap_or(b".");
        writer.write_all(read)?;
        writer.write_all(b"\t")?;
        for f in rest {
            writer.write_all(f)?;
            writer.write_all(b"\t")?;
        }
        writer.write_all(cb)?;
        writer.write_all(b"\n")?;
    }
    writer.flush()?;
    fs::rename(&tmp, pairs)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn align_chunks(
    c: &DropletConfig,
    index_dir: &Path,
    sidecar_dir: &Path,
    chunks: usize,
    join_cb: bool,
    keep_bam: bool,
    chunk_dir: &Path,
    logs: &Path,
) -> Result<Vec<ChunkResult>> {
    let next = AtomicUsize::new(0);
    let results = Mutex::new(Vec::with_capacity(chunks));
    let errors = Mutex::new(Vec::new());
    thread::scope(|scope| {
        for _ in 0..c.jobs.min(chunks) {
            let next = &next;
            let results = &results;
            let errors = &errors;
            scope.spawn(move || loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= chunks {
                    break;
                }
                match align_chunk(
                    c,
                    index_dir,
                    sidecar_dir,
                    i,
                    join_cb,
                    keep_bam,
                    chunk_dir,
                    logs,
                ) {
                    Ok(r) => results.lock().unwrap().push((i, r)),
                    Err(e) => errors.lock().unwrap().push(format!("chunk {i}: {e:#}")),
                }
            });
        }
    });
    let errors = errors.into_inner().unwrap();
    if !errors.is_empty() {
        bail!("{}", errors.join("; "))
    }
    let mut results = results.into_inner().unwrap();
    results.sort_by_key(|(i, _)| *i);
    Ok(results.into_iter().map(|(_, r)| r).collect())
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn align_chunk(
    c: &DropletConfig,
    index_dir: &Path,
    sidecar_dir: &Path,
    chunk: usize,
    join_cb: bool,
    keep_bam: bool,
    chunk_dir: &Path,
    logs: &Path,
) -> Result<ChunkResult> {
    let started = Instant::now();
    let dir = chunk_dir.join(format!("chunk-{chunk:04}"));
    fs::create_dir_all(&dir)?;
    let map = if join_cb {
        load_barcodes(&sidecar_dir.join(format!("chunk-{chunk:04}.tsv")))
            .with_context(|| format!("load barcode sidecar for chunk {chunk}"))?
    } else {
        HashMap::new()
    };
    eprintln!("event=chunk_start chunk={chunk} barcodes={}", map.len());
    let pairs_path = dir.join("pairs");
    let bam_path = dir.join("alignments.bam");
    let cancel = std::sync::Arc::new(AtomicBool::new(false));

    let fastp_r1 = dir.join("fastp_read1.fifo");
    let fastp_r3 = dir.join("fastp_read3.fifo");
    for path in [&fastp_r1, &fastp_r3] {
        let _ = fs::remove_file(path);
        fifo(path)?;
    }
    let mut fastp_command = Command::new(&c.fastp_executable);
    fastp_command
        .arg("-i")
        .arg(&fastp_r1)
        .arg("-I")
        .arg(&fastp_r3)
        .arg("--stdout")
        .arg("--adapter_sequence")
        .arg(&c.adapter)
        .arg("--adapter_sequence_r2")
        .arg(&c.adapter)
        .arg("-w")
        .arg(c.fastp_threads.to_string())
        .arg("--json")
        .arg(logs.join(format!("fastp-{chunk:04}.json")))
        .arg("--html")
        .arg(logs.join(format!("fastp-{chunk:04}.html")));
    let mut fastp = fastp_command
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::from(File::create(
            logs.join(format!("fastp-{chunk:04}.log")),
        )?))
        .spawn()
        .context("start fastp")?;
    let fastp_stdout = fastp.stdout.take().context("capture fastp stdout")?;

    let mut bwa_command = Command::new(&c.aligner_executable);
    bwa_command.arg("mem");
    if c.aligner == "bwa" {
        bwa_command.args(["-5SP", "-M", "-p"]);
    } else {
        bwa_command.args(["--hic", "-M", "-p"]);
    }
    let mut bwa = bwa_command
        .args(["-t", &c.threads_per_job.to_string()])
        .arg(&c.aligner_index)
        .arg("-")
        .stdin(Stdio::from(fastp_stdout))
        .stdout(Stdio::piped())
        .stderr(Stdio::from(File::create(
            logs.join(format!("bwa-{chunk:04}.log")),
        )?))
        .spawn()
        .context("start aligner")?;
    let sam_output = bwa.stdout.take().context("capture aligner SAM")?;

    let (bam_reader, bam_writer) = UnixStream::pair()?;
    let bam_path_thread = bam_path.clone();
    let bam_thread = if keep_bam {
        Some(thread::spawn(move || -> Result<()> {
            bam::sam_to_bam(
                BufReader::with_capacity(1024 * 1024, bam_reader),
                &bam_path_thread,
            )
        }))
    } else {
        drop(bam_reader);
        None
    };

    let (pairs_reader, pairs_writer) = UnixStream::pair()?;
    let pairs_path_thread = pairs_path.clone();
    let (min_mapq, max_mol, max_gap) = (
        c.sam.min_mapq,
        c.sam.max_molecule_size,
        c.sam.max_inter_align_gap,
    );
    let walks = c.walks_policy.clone();
    let pairs_thread = thread::spawn(move || -> Result<sam::Counts> {
        let reader = BufReader::with_capacity(1024 * 1024, pairs_reader);
        let writer = BufWriter::with_capacity(1024 * 1024, File::create(&pairs_path_thread)?);
        if walks == "all" {
            sam::convert_all_cb(reader, writer, min_mapq, max_mol, max_gap)
        } else {
            sam::convert_cb(reader, writer, min_mapq, max_mol, max_gap)
        }
    });

    let cancel_fastq_reader = std::sync::Arc::clone(&cancel);
    let index_dir_thread = index_dir.to_path_buf();
    let reader_threads = c.fastq_reader_threads;
    let fastq_reader = thread::spawn(move || -> Result<[u64; 2]> {
        stream_split(
            &index_dir_thread,
            chunk,
            [&fastp_r1, &fastp_r3],
            &cancel_fastq_reader,
            reader_threads,
        )
    });

    let mut stamper_result = stamp_sam(
        sam_output,
        &map,
        if keep_bam { Some(bam_writer) } else { None },
        pairs_writer,
        &cancel,
    );
    if stamper_result.is_err() {
        cancel.store(true, Ordering::Release);
        let _ = fastp.kill();
        let _ = bwa.kill();
    }
    let fastq_reader_result = fastq_reader
        .join()
        .unwrap_or_else(|_| Err(anyhow::anyhow!("FASTQ fastq_reader thread panicked")));
    if let Some(t) = bam_thread {
        if let Ok(Err(e)) = t.join() {
            if stamper_result.is_ok() {
                stamper_result = Err(e).context("BAM writer");
            }
        }
    }
    let pairs_result = pairs_thread
        .join()
        .unwrap_or_else(|_| Err(anyhow::anyhow!("pairs parser thread panicked")));
    let fastp_status = fastp.wait().context("wait for fastp")?;
    let bwa_status = bwa.wait().context("wait for bwa-mem3")?;
    fastq_reader_result.context("stream R1/R3")?;
    if !fastp_status.success() {
        bail!("fastp failed: {fastp_status}")
    }
    if !bwa_status.success() {
        bail!("bwa-mem3 failed: {bwa_status}")
    }
    let counts = pairs_result?;
    stamper_result?;
    let barcode_missing = missing_count();
    let result = ChunkResult {
        pairs: pairs_path,
        bam: keep_bam.then_some(bam_path),
        input_pairs: counts.input_pairs,
        selected_pairs: counts.selected_pairs,
        barcode_missing,
    };
    eprintln!(
        "event=chunk_done chunk={chunk} input_pairs={} selected_pairs={} elapsed_seconds={:.1}",
        result.input_pairs,
        result.selected_pairs,
        started.elapsed().as_secs_f64()
    );
    Ok(result)
}

thread_local! {
    static MISSING: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

pub(crate) fn missing_count() -> u64 {
    MISSING.with(|m| m.get())
}

fn stamp_sam(
    input: std::process::ChildStdout,
    map: &HashMap<Vec<u8>, Vec<u8>>,
    mut bam: Option<UnixStream>,
    mut pairs: UnixStream,
    cancel: &AtomicBool,
) -> Result<()> {
    let mut reader = BufReader::with_capacity(1024 * 1024, input);
    let mut line = Vec::new();
    let mut out = Vec::with_capacity(4096);
    let mut body = false;
    let mut missing = 0u64;
    let mut lines = 0u64;
    loop {
        if cancel.load(Ordering::Relaxed) {
            bail!("cancelled")
        }
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if !body && line.starts_with(b"@") {
            if let Some(w) = bam.as_mut() {
                w.write_all(&line)?;
            }
            pairs.write_all(&line)?;
            continue;
        }
        body = true;
        let end = line.len() - usize::from(line.last() == Some(&b'\n'));
        let end = if line.get(end.wrapping_sub(1)) == Some(&b'\r') {
            end - 1
        } else {
            end
        };
        let qname_end = line[..end].iter().position(|b| *b == b'\t').unwrap_or(end);
        let cb = map.get(&line[..qname_end]);
        if cb.is_none() {
            missing += 1;
        }
        out.clear();
        out.extend_from_slice(&line[..end]);
        out.extend_from_slice(b"\tCB:Z:");
        out.extend_from_slice(cb.map(Vec::as_slice).unwrap_or(b"."));
        out.push(b'\n');
        if let Some(w) = bam.as_mut() {
            w.write_all(&out)?;
        }
        pairs.write_all(&out)?;
        lines += 1;
        if lines.is_multiple_of(4096) {
            if let Some(w) = bam.as_mut() {
                w.flush()?;
            }
            pairs.flush()?;
        }
    }
    if let Some(mut w) = bam {
        w.flush()?;
    }
    pairs.flush()?;
    MISSING.with(|m| m.set(missing));
    if missing > 0 {
        eprintln!("warning: {missing} SAM records had no R2 barcode entry");
    }
    Ok(())
}

fn merge_pairs(inputs: &[PathBuf], output: &Path) -> Result<()> {
    let mut out = BufWriter::with_capacity(1024 * 1024, File::create(output)?);
    out.write_all(PAIRS_HEADER.as_bytes())?;
    for path in inputs {
        for line in BufReader::with_capacity(1024 * 1024, File::open(path)?).lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let kind = line
                .split('\t')
                .nth(7)
                .with_context(|| format!("pairs row has no pair_type in {}", path.display()))?;
            if !matches!(kind, "UU" | "UR" | "RU") {
                continue;
            }
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

fn sort_pairs(c: &DropletConfig, input: &Path, output: &Path, logs: &Path) -> Result<()> {
    let status = stage_command(Stage::SortPairs)?
        .arg(&c.chromsizes)
        .arg(&c.assembly)
        .arg(output)
        .arg(c.scratch())
        .arg(c.sort_threads.to_string())
        .arg(c.sort_run_mib.to_string())
        .arg(input)
        .stderr(Stdio::from(File::create(logs.join("sort.log"))?))
        .status()
        .context("run pair sort")?;
    if !status.success() {
        bail!(
            "pair sort failed: {status}; inspect {}",
            logs.join("sort.log").display()
        )
    }
    Ok(())
}

fn build_matrix(c: &DropletConfig, pairs: &Path, logs: &Path) -> Result<()> {
    let matrix_dir = c.output_dir.join("matrix");
    fs::create_dir_all(&matrix_dir)?;
    let base = c.scratch().join("base.cool");
    if base.exists() {
        fs::remove_file(&base)?;
    }
    let mut bin = stage_command(Stage::BinPairs)?;
    bin.arg(&c.chromsizes)
        .arg(c.matrix.resolutions[0].to_string())
        .arg(&c.assembly)
        .arg(pairs)
        .arg(&base)
        .arg(c.scratch())
        .arg(c.sort_threads.to_string());
    crate::matrix::chromosome_options(&mut bin, c.matrix.sort_chromosomes);
    let status = bin
        .stderr(Stdio::from(File::create(logs.join("matrix-bin.log"))?))
        .status()
        .context("run pairs-to-cool")?;
    if !status.success() {
        bail!(
            "pairs-to-cool failed: {status}; inspect {}",
            logs.join("matrix-bin.log").display()
        )
    }
    let mcool = matrix_dir.join(format!("{}.mcool", c.prefix));
    if mcool.exists() {
        fs::remove_file(&mcool)?;
    }
    let resolutions = c
        .matrix
        .resolutions
        .iter()
        .map(u64::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let mut multi = stage_command(Stage::BuildMcool)?;
    multi
        .arg(&base)
        .arg(&c.chromsizes)
        .arg(&mcool)
        .arg(&c.assembly)
        .arg(resolutions)
        .arg(c.sort_threads.to_string())
        .env(
            crate::matrix::BALANCE_MAX_ITERS_ENV,
            c.matrix.balance_max_iterations.to_string(),
        );
    crate::matrix::chromosome_options(&mut multi, c.matrix.sort_chromosomes);
    let status = multi
        .stderr(Stdio::from(File::create(logs.join("matrix-mcool.log"))?))
        .status()
        .context("run cool-to-mcool")?;
    if !status.success() {
        bail!(
            "cool-to-mcool failed: {status}; inspect {}",
            logs.join("matrix-mcool.log").display()
        )
    }
    eprintln!("event=matrix_done path={}", mcool.display());
    Ok(())
}

fn dedup_pairs(c: &DropletConfig, input: &Path, output: &Path, logs: &Path) -> Result<()> {
    let status = stage_command(Stage::DedupPairs)?
        .arg(c.dedup.workers.to_string())
        .arg(c.dedup.max_mismatch.to_string())
        .stdin(Stdio::from(File::open(input)?))
        .stdout(Stdio::from(File::create(output)?))
        .stderr(Stdio::from(File::create(logs.join("dedup.log"))?))
        .status()
        .context("run pair dedup")?;
    if !status.success() {
        bail!(
            "pair dedup failed: {status}; inspect {}",
            logs.join("dedup.log").display()
        )
    }
    Ok(())
}

fn write_pipeline_qc(
    qc_dir: &Path,
    summary: &qc::Summary,
    input_pairs: u64,
    selected_pairs: u64,
    barcode_missing: u64,
    seconds: f64,
) -> Result<()> {
    #[derive(serde::Serialize)]
    struct PipelineQc<'a> {
        input_pairs: u64,
        selected_pairs: u64,
        barcode_missing_sam_records: u64,
        elapsed_seconds: f64,
        qc: &'a qc::Summary,
    }
    fs::write(
        qc_dir.join("pipeline.json"),
        serde_json::to_vec_pretty(&PipelineQc {
            input_pairs,
            selected_pairs,
            barcode_missing_sam_records: barcode_missing,
            elapsed_seconds: seconds,
            qc: summary,
        })?,
    )?;
    Ok(())
}

fn write_manifest(
    c: &DropletConfig,
    pairs: Option<&Path>,
    pairs_gz: Option<&Path>,
    summary: &qc::Summary,
) -> Result<()> {
    let mut out = BufWriter::new(File::create(c.output_dir.join("MANIFEST.tsv"))?);
    writeln!(out, "kind\tpath")?;
    if let Some(path) = pairs {
        writeln!(out, "pairs\t{}", path.display())?;
    }
    if let Some(path) = pairs_gz {
        writeln!(out, "pairs_gz\t{}", path.display())?;
    }
    if c.keep_bam {
        writeln!(
            out,
            "bam\t{}",
            c.output_dir
                .join("alignments")
                .join(format!("{}.bam", c.prefix))
                .display()
        )?;
    }
    writeln!(
        out,
        "qc_per_cell\t{}",
        c.output_dir.join("qc/per_cell.tsv").display()
    )?;
    if c.matrix.enabled {
        writeln!(
            out,
            "mcool\t{}",
            c.output_dir
                .join("matrix")
                .join(format!("{}.mcool", c.prefix))
                .display()
        )?;
    }
    if c.atac_fragments {
        writeln!(
            out,
            "fragments\t{}",
            c.output_dir
                .join("atac")
                .join(format!("{}.fragments.tsv.gz", c.prefix))
                .display()
        )?;
    }
    writeln!(
        out,
        "qc_barcodes\t{}",
        c.output_dir.join("qc/barcode_stats.tsv").display()
    )?;
    writeln!(out, "pairs_total\t{}", summary.total_pairs)?;
    writeln!(out, "cells\t{}", summary.cells)?;
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn join_pairs_replaces_cb_and_keeps_nine_columns() {
        let dir = std::env::temp_dir().join(format!("rh_join_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let pairs = dir.join("pairs");
        std::fs::write(
            &pairs,
            "## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type CB\nr1\tchr1\t10\tchr1\t20\t+\t-\tUU\t.\nr2\tchr1\t30\tchr2\t40\t+\t-\tUR\t.\n",
        )
        .unwrap();
        let sidecar = dir.join("chunk-0000.tsv");
        std::fs::write(&sidecar, "r1\tAAAACCCCGGGGTTTT\n").unwrap();
        join_pairs_cb(&pairs, &sidecar).unwrap();
        let text = std::fs::read_to_string(&pairs).unwrap();
        let body: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(body.len(), 2);
        let f: Vec<&str> = body[0].split('\t').collect();
        assert_eq!(f.len(), 9);
        assert_eq!(f[0], "r1");
        assert_eq!(f[1], "chr1");
        assert_eq!(f[8], "AAAACCCCGGGGTTTT");
        let g: Vec<&str> = body[1].split('\t').collect();
        assert_eq!(g.len(), 9);
        assert_eq!(g[8], ".");
        std::fs::remove_dir_all(&dir).ok();
    }
}
