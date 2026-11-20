//! Slurm job-array execution for the droplet (3-read) pipeline.
//!
//! Phases: `prepare` (FASTQ index + barcode prepass), `chunk` (one array task
//! per chunk: fastp + align + parse), and `finish` (merge, sort, per-cell
//! dedup, QC, matrix). The scheduler plumbing — `sbatch`, scripts, logs,
//! `slurm-jobs.json` — is the same code the bulk v0.5.0 pipeline uses.
use super::{barcode, config::DropletConfig, index::IndexManifest, run};
use crate::{
    cli::run::SlurmOptions,
    workflow::slurm::{
        cancel_submitted, clear_failure, now_ms, record_failure, release_array,
        sbatch_named, validate_extra_args, write_json, write_task_script, Jobs, SCHEMA,
    },
};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    ffi::OsStr,
    fs::{self, File},
    os::unix::fs::PermissionsExt,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

#[derive(Serialize, Deserialize)]
pub(crate) struct SlurmPlan {
    schema: u32,
    run_id: String,
    chunks: usize,
    keep_bam: bool,
    output: PathBuf,
    work: PathBuf,
    config: DropletConfig,
}

#[derive(Serialize, Deserialize)]
struct ChunkMarker {
    schema: u32,
    chunk: usize,
    start_ms: u64,
    end_ms: u64,
    input_pairs: u64,
    selected_pairs: u64,
    barcode_missing: u64,
    pairs_bytes: u64,
    bam_bytes: Option<u64>,
}

/// The droplet indexer derives its chunk count from compressed sizes, so the
/// array size is fixed at submission time. An existing index is reused.
fn plan_chunks(c: &DropletConfig) -> Result<usize> {
    let dir = c.index();
    if dir.join("manifest.json").is_file() {
        return Ok(IndexManifest::load(&dir)?.chunks);
    }
    let size = crate::cli::index::parse_chunk_size(OsStr::new(&c.chunk_size))?;
    let r1 = fs::metadata(&c.read1)?.len();
    let mate2 = fs::metadata(c.mate2()?)?.len();
    let chunks = usize::try_from(r1.max(mate2).div_ceil(size))?;
    if !(1..=256).contains(&chunks) {
        bail!("estimated droplet chunk count {chunks} is outside 1..=256; increase chunk_size")
    }
    Ok(chunks)
}

fn job_cpus(c: &DropletConfig, args: &SlurmOptions) -> (usize, usize, usize) {
    let prepare = args
        .index_cpus
        .unwrap_or_else(|| 8usize.max(c.barcode.threads as usize))
        .max(1);
    let chunk = args
        .map_cpus
        .unwrap_or(c.threads_per_job + c.fastp_threads + 4)
        .max(1);
    let finish_base = c
        .sort_threads
        .max(c.dedup.workers)
        .max(if c.matrix.enabled { 4 } else { 1 });
    let finish = args.finish_cpus.unwrap_or(finish_base + 4).max(1);
    (prepare, chunk, finish)
}

fn ensure_sbatch() -> Result<()> {
    let status = std::process::Command::new("sbatch")
        .arg("--version")
        .output()
        .context("start sbatch; is Slurm available on PATH?")?;
    if !status.status.success() {
        bail!("sbatch --version failed")
    }
    Ok(())
}

pub(crate) fn submit(config: &Path, args: SlurmOptions, check: bool, keep_bam: bool) -> Result<()> {
    let (mut c, _base) = DropletConfig::load(config)?;
    c.keep_bam |= keep_bam;
    let chunks = plan_chunks(&c)?;
    let parallel = args.max_parallel.unwrap_or(c.jobs).clamp(1, chunks);
    let (prepare_cpus, chunk_cpus, finish_cpus) = job_cpus(&c, &args);
    let time = args.time.as_deref().unwrap_or("04:00:00");
    let prepare_mem = args.index_mem.as_deref().unwrap_or("16G");
    let chunk_mem = args.map_mem.as_deref().unwrap_or("32G");
    let finish_mem = args.finish_mem.as_deref().unwrap_or("64G");

    if check {
        println!(
            "Droplet Slurm configuration valid\n  Output: {}\n  Chunks: {} (up to {} concurrent)\n  CPUs per job: prepare {}, chunk {}, finish {}\n  Memory per job: prepare {}, chunk {}, finish {}\n  Phases: prepare (index + {}) -> chunk array ({}) -> finish (merge/sort/dedup/QC/matrix)",
            c.output_dir.display(),
            chunks,
            parallel,
            prepare_cpus,
            chunk_cpus,
            finish_cpus,
            prepare_mem,
            chunk_mem,
            finish_mem,
            if c.is_two_read() {
                "barcode orientation"
            } else {
                "barcode prepass"
            },
            if c.is_two_read() {
                "fused barcode + align"
            } else {
                "align"
            }
        );
        return Ok(());
    }

    validate_extra_args(&args.sbatch_arg)?;
    ensure_sbatch()?;

    let out = c.output_dir.clone();
    fs::create_dir_all(&out).with_context(|| format!("create {}", out.display()))?;
    let work = out.join("slurm");
    fs::create_dir_all(work.join("chunks"))?;
    let jobs_file = out.join("slurm-jobs.json");
    if jobs_file.exists() {
        bail!(
            "Slurm job record already exists: {}; choose a new output_dir",
            jobs_file.display()
        )
    }
    let logs = out.join("logs").join("slurm");
    fs::create_dir_all(&logs)?;
    let exe = work.join("rust-hic");
    fs::copy(std::env::current_exe()?, &exe).context("copy executable to shared work directory")?;
    fs::set_permissions(&exe, fs::Permissions::from_mode(0o755))?;

    let run_id = format!("{:08x}", now_ms()? as u32);
    let plan_path = work.join("plan.yaml");
    let plan = SlurmPlan {
        schema: SCHEMA,
        run_id: run_id.clone(),
        chunks,
        keep_bam: c.keep_bam,
        output: out.clone(),
        work: work.clone(),
        config: c,
    };
    serde_yaml::to_writer(File::create(&plan_path)?, &plan)?;
    fs::write(out.join("work-dir.txt"), format!("{}\n", work.display()))?;

    write_task_script(
        &work.join("prepare.sh"),
        &exe,
        "droplet-slurm-task",
        "prepare",
        &plan_path,
        false,
    )?;
    write_task_script(
        &work.join("chunk.sh"),
        &exe,
        "droplet-slurm-task",
        "chunk",
        &plan_path,
        true,
    )?;
    write_task_script(
        &work.join("finish.sh"),
        &exe,
        "droplet-slurm-task",
        "finish",
        &plan_path,
        false,
    )?;

    let group_id = run_id.clone();
    let mut jobs = Jobs {
        schema: SCHEMA,
        run_id,
        group_id: Some(group_id),
        chunks,
        index: None,
        alignment: None,
        finish: None,
    };
    write_json(&jobs_file, &jobs)?;

    let mut submitted = Vec::new();
    let result = (|| -> Result<()> {
        let id = sbatch_named(
            &args,
            "rust-hic-droplet-prepare",
            "droplet prepare",
            &logs.join("prepare-%j.log"),
            prepare_cpus,
            prepare_mem,
            args.index_time.as_deref().unwrap_or(time),
            None,
            None,
            &work.join("prepare.sh"),
            false,
        )?;
        submitted.push(id.clone());
        jobs.index = Some(id.clone());
        write_json(&jobs_file, &jobs)?;

        let array = format!("0-{}%{parallel}", chunks - 1);
        let id = sbatch_named(
            &args,
            "rust-hic-droplet-chunk",
            "droplet chunk array",
            &logs.join("chunk-%A_%a.log"),
            chunk_cpus,
            chunk_mem,
            args.map_time.as_deref().unwrap_or(time),
            Some(&format!("afterok:{}", jobs.index.as_deref().unwrap_or_default())),
            Some(&array),
            &work.join("chunk.sh"),
            false,
        )?;
        submitted.push(id.clone());
        jobs.alignment = Some(id.clone());
        write_json(&jobs_file, &jobs)?;
        // The shared submitter holds arrays until task names are assigned.
        release_array(&id)?;

        let id = sbatch_named(
            &args,
            "rust-hic-droplet-finish",
            "droplet final processing",
            &logs.join("finish-%j.log"),
            finish_cpus,
            finish_mem,
            args.finish_time.as_deref().unwrap_or(time),
            Some(&format!("afterany:{id}")),
            None,
            &work.join("finish.sh"),
            true,
        )?;
        submitted.push(id.clone());
        jobs.finish = Some(id);
        write_json(&jobs_file, &jobs)?;
        Ok(())
    })();
    if let Err(e) = result {
        if !submitted.is_empty() {
            let _ = cancel_submitted(&submitted);
        }
        record_failure(&out, "submission", &e);
        return Err(e).with_context(|| format!("Slurm submission failed; inspect {}", out.display()));
    }

    println!(
        "Submitted droplet Slurm workflow\n  Output: {}\n  Prepare job: {}\n  Chunk array: {} ({} chunks, up to {} concurrent)\n  Finish job: {}\n  Status: rust-hic run status {}\n  Logs: {}",
        out.display(),
        jobs.index.as_deref().unwrap_or("?"),
        jobs.alignment.as_deref().unwrap_or("?"),
        chunks,
        parallel,
        jobs.finish.as_deref().unwrap_or("?"),
        out.display(),
        logs.display()
    );
    Ok(())
}

/// Worker entry point invoked inside the Slurm allocation.
pub(crate) fn task(phase: &str, plan_path: &Path, chunk: Option<usize>) -> Result<()> {
    let plan: SlurmPlan = serde_yaml::from_reader(
        File::open(plan_path).with_context(|| format!("open {}", plan_path.display()))?,
    )
    .with_context(|| format!("parse {}", plan_path.display()))?;
    if plan.schema != SCHEMA {
        bail!("unsupported droplet Slurm plan schema {}", plan.schema)
    }
    clear_failure(&plan.output, phase)?;
    let result = match phase {
        "prepare" => prepare(&plan),
        "chunk" => chunk_task(
            &plan,
            chunk.context("droplet Slurm chunk task requires a chunk index")?,
        ),
        "finish" => finish(&plan),
        other => bail!("unknown droplet Slurm phase: {other}"),
    };
    if let Err(error) = &result {
        record_failure(&plan.output, phase, error);
    }
    result
}

fn prepare(plan: &SlurmPlan) -> Result<()> {
    let c = &plan.config;
    run::create_dirs(c)?;
    let logs = c.output_dir.join("logs");
    let index_dir = crate::droplet::index::ensure_index(c, &logs.join("fastq-index.log"))?;
    let manifest = IndexManifest::load(&index_dir)?;
    if manifest.chunks != plan.chunks {
        bail!(
            "FASTQ index produced {} chunks but the Slurm array was submitted with {}; resubmit to rebuild the array",
            manifest.chunks,
            plan.chunks
        )
    }
    if c.is_two_read() {
        let chosen = crate::droplet::two_read::resolve_orientation(c, &index_dir)?;
        eprintln!(
            "event=droplet_slurm_prepare chunks={} barcode_orientation={chosen}",
            manifest.chunks
        );
    } else {
        let prepass = barcode::prepass(c, &manifest, &c.scratch().join("barcodes"))?;
        run::write_barcode_qc(&prepass, &c.output_dir.join("qc"))?;
        eprintln!(
            "event=droplet_slurm_prepare chunks={} assignable={}",
            manifest.chunks,
            prepass
                .total
                .exact
                .saturating_add(prepass.total.corrected)
        );
    }
    Ok(())
}

fn chunk_task(plan: &SlurmPlan, chunk: usize) -> Result<()> {
    let c = &plan.config;
    let logs = c.output_dir.join("logs");
    let index_dir = c.index();
    let manifest = IndexManifest::load(&index_dir)?;
    if chunk >= manifest.chunks {
        bail!(
            "chunk {chunk} is outside the index ({} chunks)",
            manifest.chunks
        )
    }
    let chunk_dir = c.scratch().join("chunks");
    fs::create_dir_all(&chunk_dir)?;
    let sidecar_dir = c.scratch().join("barcodes");
    let start_ms = now_ms()?;
    let result = run::align_chunk(
        c,
        &index_dir,
        &sidecar_dir,
        chunk,
        !c.is_two_read() && c.cb_join != "post",
        plan.keep_bam,
        &chunk_dir,
        &logs,
    )?;
    let end_ms = now_ms()?;
    let marker = ChunkMarker {
        schema: SCHEMA,
        chunk,
        start_ms,
        end_ms,
        input_pairs: result.input_pairs,
        selected_pairs: result.selected_pairs,
        barcode_missing: result.barcode_missing,
        pairs_bytes: fs::metadata(&result.pairs)?.len(),
        bam_bytes: result
            .bam
            .as_ref()
            .and_then(|p| fs::metadata(p).ok())
            .map(|m| m.len()),
    };
    let dir = plan.work.join("chunks").join(format!("{chunk:04}"));
    fs::create_dir_all(&dir)?;
    let mut file = File::create(dir.join("complete.json"))?;
    serde_json::to_writer_pretty(&mut file, &marker)?;
    Ok(())
}

fn finish(plan: &SlurmPlan) -> Result<()> {
    let c = &plan.config;
    let chunk_dir = c.scratch().join("chunks");
    let mut results = Vec::with_capacity(plan.chunks);
    let mut missing = Vec::new();
    let (mut earliest, mut latest) = (u64::MAX, 0u64);
    for chunk in 0..plan.chunks {
        let path = plan
            .work
            .join("chunks")
            .join(format!("{chunk:04}"))
            .join("complete.json");
        let marker: ChunkMarker = match File::open(&path) {
            Ok(file) => {
                serde_json::from_reader(file).with_context(|| format!("parse {}", path.display()))?
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                missing.push(chunk);
                continue;
            }
            Err(error) => return Err(error).with_context(|| format!("open {}", path.display())),
        };
        earliest = earliest.min(marker.start_ms);
        latest = latest.max(marker.end_ms);
        let bam = if plan.keep_bam {
            let path = chunk_dir
                .join(format!("chunk-{chunk:04}"))
                .join("alignments.bam");
            if !path.is_file() {
                bail!(
                    "chunk {chunk}: keep_bam is enabled but {} is missing",
                    path.display()
                )
            }
            Some(path)
        } else {
            None
        };
        results.push(run::ChunkResult {
            pairs: chunk_dir.join(format!("chunk-{chunk:04}")).join("pairs"),
            bam,
            input_pairs: marker.input_pairs,
            selected_pairs: marker.selected_pairs,
            barcode_missing: marker.barcode_missing,
        });
    }
    if !missing.is_empty() {
        bail!("missing chunk results for {missing:?}; refusing to finish")
    }
    let total = Duration::from_millis(latest.saturating_sub(earliest));
    let started = Instant::now()
        .checked_sub(total)
        .unwrap_or_else(Instant::now);
    run::finish_pipeline(c, started, &results, plan.keep_bam)
}
