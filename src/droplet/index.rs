//! Reuse the v1 paired FASTQ index for R1/R3 and read its boundaries.
use anyhow::{bail, Context, Result};
use serde::Deserialize;
use std::{
    fs::{self, File},
    path::{Path, PathBuf},
    process::Command,
};

use super::config::DropletConfig;

#[derive(Debug, Deserialize)]
pub(crate) struct Boundary {
    pub read_id: String,
}

#[derive(Debug, Deserialize)]
pub(crate) struct IndexManifest {
    pub version: u32,
    pub chunks: usize,
    pub boundaries: Vec<Boundary>,
}

impl IndexManifest {
    pub(crate) fn load(dir: &Path) -> Result<Self> {
        let path = dir.join("manifest.json");
        let m: Self = serde_json::from_reader(
            File::open(&path).with_context(|| format!("open {}", path.display()))?,
        )
        .with_context(|| format!("parse {}", path.display()))?;
        if m.version != 1 || m.chunks == 0 || m.boundaries.len() != m.chunks + 1 {
            bail!("unsupported FASTQ index manifest {}", path.display())
        }
        Ok(m)
    }
    /// Read IDs that start each chunk after the first.
    pub(crate) fn chunk_starts(&self) -> Vec<&str> {
        self.boundaries[1..self.chunks]
            .iter()
            .map(|b| b.read_id.as_str())
            .collect()
    }
}

/// Build or reuse the R1/R3 index. The index stores raw input fingerprints;
/// barcode work never mutates the FASTQs.
pub(crate) fn ensure_index(c: &DropletConfig, log: &Path) -> Result<PathBuf> {
    let dir = c.index();
    if dir.join("manifest.json").is_file() {
        eprintln!("event=fastq_index_reuse path={}", dir.display());
        return Ok(dir);
    }
    if dir.exists() {
        bail!(
            "FASTQ index directory exists without manifest.json: {}; remove it or point index_dir elsewhere",
            dir.display()
        )
    }
    if let Some(parent) = dir.parent() {
        fs::create_dir_all(parent)?;
    }
    eprintln!("event=fastq_index_start path={}", dir.display());
    let exe = std::env::current_exe().context("resolve rust-hic executable")?;
    let status = Command::new(exe)
        .arg("index")
        .arg(&c.read1)
        .arg(c.mate2()?)
        .arg("-o")
        .arg(&dir)
        .arg("--chunk-size")
        .arg(&c.chunk_size)
        .arg("--threads-per-mate")
        .arg("2")
        .arg("--search-mib")
        .arg(c.index_search_mib.to_string())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::from(File::create(log)?))
        .status()
        .context("run FASTQ indexer")?;
    if !status.success() {
        bail!("FASTQ indexer failed: {status}; inspect {}", log.display())
    }
    let manifest = IndexManifest::load(&dir)?;
    eprintln!(
        "event=fastq_index_done chunks={} path={}",
        manifest.chunks,
        dir.display()
    );
    Ok(dir)
}
