//! Configuration for the droplet (3-read) Hi-C pipeline.
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    fs::{self, File},
    path::{Path, PathBuf},
};

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct DropletConfig {
    /// R1: gDNA (Nextera read 1 side).
    pub read1: PathBuf,
    /// R2: cell barcode read.
    pub barcode_read: PathBuf,
    /// R3: gDNA mate.
    pub read3: PathBuf,
    /// Reference FASTA (same path as the aligner index prefix).
    pub reference: PathBuf,
    /// bwa-mem3/bwa-mem2 index prefix.
    pub aligner_index: PathBuf,
    /// Aligner executable (bwa-mem3 recommended).
    #[serde(default = "default_bwa_mem3")]
    pub aligner_executable: PathBuf,
    /// Aligner flag style: `bwa-mem3` (`--hic`) or `bwa` (`-5SP`).
    #[serde(default = "default_aligner")]
    pub aligner: String,
    /// Chromosome sizes for pair sorting and QC.
    pub chromsizes: PathBuf,
    /// 16 bp barcode whitelist (one barcode per line).
    pub whitelist: PathBuf,
    /// New or reusable output directory.
    pub output_dir: PathBuf,
    /// Output prefix (matrix/pairs/QC names).
    pub prefix: String,
    /// fastp executable.
    #[serde(default = "default_fastp")]
    pub fastp_executable: PathBuf,
    /// Adapter sequence for both gDNA reads (Nextera mosaic end).
    #[serde(default = "default_adapter")]
    pub adapter: String,
    /// Scratch directory; defaults to <output_dir>/scratch.
    #[serde(default)]
    pub scratch_dir: Option<PathBuf>,
    /// Reusable FASTQ index directory; defaults to <output_dir>/fastq-index.
    #[serde(default)]
    pub index_dir: Option<PathBuf>,
    /// Number of chunks processed concurrently.
    #[serde(default = "one")]
    pub jobs: usize,
    /// Alignment threads per chunk (bwa-mem3 -t).
    #[serde(default = "eight")]
    pub threads_per_job: usize,
    /// fastp worker threads per chunk.
    #[serde(default = "two")]
    pub fastp_threads: usize,
    /// Decompression threads for reading the R1/R3 FASTQ streams.
    #[serde(default = "default_fastq_reader_threads", alias = "fastq_reader_threads")]
    pub fastq_reader_threads: u32,
    /// Target compressed chunk size accepted by the FASTQ indexer.
    #[serde(default = "default_chunk_size")]
    pub chunk_size: String,
    /// Boundary search window in MiB (increase for concatenated runs).
    #[serde(default = "default_search_mib")]
    pub index_search_mib: usize,
    /// Retain BAM (per-chunk BAMs are removed after a successful merge).
    #[serde(default = "no")]
    pub keep_bam: bool,
    /// Retain the final pairs as a deliverable. The format is plain text unless
    /// `compressed_pairs` is set.
    #[serde(default = "yes")]
    pub keep_pairs: bool,
    /// Retain the pairs as BGZF `<prefix>.pairs.gz` with a pairix `.px2` index.
    /// This supersedes the plain file: with `keep_pairs: true` and
    /// `compressed_pairs: true` only the compressed file is kept.
    #[serde(default)]
    pub compressed_pairs: bool,
    /// How cell barcodes are attached: `fused` stamps them during parsing,
    /// `post` aligns first and joins afterwards, `none` has no barcode read.
    /// `cb_join` remains accepted as an alias for existing configs.
    #[serde(default = "default_cell_barcode_process_mode", alias = "cb_join")]
    pub cell_barcode_process_mode: String,
    /// In `post` mode, run the barcode preprocessing concurrently with alignment.
    #[serde(default, alias = "prepass_overlap")]
    pub barcode_preprocessing_overlap: bool,
    /// Reuse barcode sidecar files from a previous run.
    #[serde(default = "yes")]
    pub reuse_barcodes: bool,
    /// Genome assembly label written to pairs headers.
    #[serde(default = "default_assembly")]
    pub assembly: String,
    /// Walk policy: `5unique` (default) or `all` (multi-way contacts).
    #[serde(default = "default_walks")]
    pub walks_policy: String,
    /// Emit Tn5-shifted ATAC fragments from the merged BAM (ATAC-Hi-C/spatial).
    #[serde(default)]
    pub atac_fragments: bool,
    #[serde(default = "four_u64")]
    pub t5_shift: u64,
    #[serde(default = "one_u64")]
    pub atac_min_length: u64,
    #[serde(default = "five_k")]
    pub atac_max_length: u64,
    /// Pair-sort threads.
    #[serde(default = "four")]
    pub sort_threads: usize,
    /// In-memory run size for pair sorting, in MiB.
    #[serde(default = "default_run_mib")]
    pub sort_run_mib: usize,
    #[serde(default)]
    pub barcode: BarcodeConfig,
    #[serde(default)]
    pub sam: SamConfig,
    #[serde(default)]
    pub dedup: DedupConfig,
    #[serde(default)]
    pub matrix: MatrixConfig,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub(crate) struct BarcodeConfig {
    /// First base of the barcode in R2 (0-based).
    pub slice_start: usize,
    /// Barcode length in bases.
    pub slice_len: usize,
    /// `rc`, `forward` or `auto`.
    pub orientation: String,
    /// Maximum Hamming distance considered (0, 1 or 2 implemented).
    pub max_mismatch: u8,
    /// Threads used to match R2 barcodes against the whitelist.
    #[serde(default = "default_barcode_threads", alias = "threads")]
    pub matcher_threads: u32,
}

impl Default for BarcodeConfig {
    fn default() -> Self {
        Self {
            slice_start: 0,
            slice_len: 16,
            orientation: "rc".into(),
            max_mismatch: 1,
            matcher_threads: default_barcode_threads(),
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub(crate) struct SamConfig {
    pub min_mapq: u8,
    pub max_molecule_size: i64,
    pub max_inter_align_gap: i64,
}

impl Default for SamConfig {
    fn default() -> Self {
        Self {
            min_mapq: 1,
            max_molecule_size: 750,
            max_inter_align_gap: 20,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub(crate) struct DedupConfig {
    pub max_mismatch: i64,
    pub workers: usize,
}

impl Default for DedupConfig {
    fn default() -> Self {
        Self {
            max_mismatch: 3,
            workers: 4,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub(crate) struct MatrixConfig {
    /// Build the pseudo-bulk multi-resolution mcool.
    pub enabled: bool,
    /// Smallest resolution first; every value must be a multiple of it.
    pub resolutions: Vec<u64>,
    pub balance_max_iterations: usize,
    pub sort_chromosomes: bool,
}

impl Default for MatrixConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            resolutions: vec![
                5_000, 10_000, 25_000, 50_000, 100_000, 250_000, 500_000, 1_000_000, 2_500_000,
                5_000_000,
            ],
            balance_max_iterations: 1_000,
            sort_chromosomes: true,
        }
    }
}

fn default_bwa_mem3() -> PathBuf {
    "/projects/b1042/YueLab/JW/rusthic/bin/bwa-mem3".into()
}
fn default_aligner() -> String {
    "bwa-mem3".into()
}
fn default_fastp() -> PathBuf {
    "/projects/b1042/YueLab/JW/rusthic-droplet/build/bin/fastp".into()
}
fn default_adapter() -> String {
    "CTGTCTCTTATACACATCT".into()
}
fn default_chunk_size() -> String {
    "2GiB".into()
}
fn default_search_mib() -> usize {
    1024
}
fn default_barcode_threads() -> u32 {
    8
}
fn default_cell_barcode_process_mode() -> String {
    "fused".into()
}
fn default_assembly() -> String {
    "hg38".into()
}
fn default_walks() -> String {
    "5unique".into()
}
fn default_run_mib() -> usize {
    256
}
fn one() -> usize {
    1
}
fn two() -> usize {
    2
}
fn default_fastq_reader_threads() -> u32 {
    4
}
fn four() -> usize {
    4
}
fn four_u64() -> u64 {
    4
}
fn one_u64() -> u64 {
    1
}
fn five_k() -> u64 {
    5_000
}
fn eight() -> usize {
    8
}
fn yes() -> bool {
    true
}
fn no() -> bool {
    false
}

fn resolve(base: &Path, p: PathBuf) -> PathBuf {
    if p.is_absolute() {
        p
    } else {
        base.join(p)
    }
}

impl DropletConfig {
    pub(crate) fn load(path: &Path) -> Result<(Self, PathBuf)> {
        let path = fs::canonicalize(path)
            .with_context(|| format!("open droplet configuration {}", path.display()))?;
        let base = path
            .parent()
            .context("configuration has no parent directory")?
            .to_path_buf();
        let mut c: Self = serde_yaml::from_reader(
            File::open(&path).with_context(|| format!("open {}", path.display()))?,
        )
        .with_context(|| format!("parse {}", path.display()))?;
        for p in [
            &mut c.read1,
            &mut c.barcode_read,
            &mut c.read3,
            &mut c.reference,
            &mut c.aligner_index,
            &mut c.aligner_executable,
            &mut c.chromsizes,
            &mut c.whitelist,
            &mut c.output_dir,
            &mut c.fastp_executable,
        ] {
            *p = resolve(&base, std::mem::take(p));
        }
        if let Some(p) = c.scratch_dir.take() {
            c.scratch_dir = Some(resolve(&base, p));
        }
        if let Some(p) = c.index_dir.take() {
            c.index_dir = Some(resolve(&base, p));
        }
        c.validate()?;
        Ok((c, base))
    }

    pub(crate) fn validate(&self) -> Result<()> {
        if self.prefix.is_empty()
            || self.prefix.contains(['/', '\\', '\0'])
            || self.prefix == "."
            || self.prefix == ".."
        {
            bail!("prefix must be a non-empty filename component")
        }
        if self.jobs == 0 || self.jobs > 64 {
            bail!("jobs must be 1..=64")
        }
        if self.threads_per_job == 0 || self.threads_per_job > 128 {
            bail!("threads_per_job must be 1..=128")
        }
        if self.fastp_threads == 0 || self.fastp_threads > 16 {
            bail!("fastp_threads must be 1..=16")
        }
        if !(1..=256).contains(&self.sort_threads) {
            bail!("sort_threads must be 1..=256")
        }
        if self.sort_run_mib == 0 {
            bail!("sort_run_mib must be positive")
        }
        if self.barcode.slice_len == 0 || self.barcode.slice_len > 16 {
            bail!("barcode.slice_len must be 1..=16")
        }
        if self.barcode.max_mismatch > 2 {
            bail!("barcode.max_mismatch currently supports 0, 1 or 2")
        }
        if self.barcode.matcher_threads == 0 || self.barcode.matcher_threads > 32 {
            bail!("barcode.matcher_threads must be 1..=32")
        }
        if self.fastq_reader_threads == 0 || self.fastq_reader_threads > 16 {
            bail!("fastq_reader_threads must be 1..=16")
        }
        match self.cell_barcode_process_mode.as_str() {
            "fused" | "post" | "none" => {}
            other => bail!(
                "cell_barcode_process_mode must be fused, post or none, got {other:?}"
            ),
        }
        match self.barcode.orientation.as_str() {
            "rc" | "forward" | "auto" => {}
            other => bail!("barcode.orientation must be rc, forward or auto, got {other:?}"),
        }
        match self.aligner.as_str() {
            "bwa-mem3" | "bwa" => {}
            other => bail!("aligner must be bwa-mem3 or bwa, got {other:?}"),
        }
        match self.walks_policy.as_str() {
            "5unique" | "all" => {}
            other => bail!("walks_policy must be 5unique or all, got {other:?}"),
        }
        if self.matrix.enabled {
            let r = &self.matrix.resolutions;
            if r.is_empty() || r[0] == 0 || r.iter().any(|x| x % r[0] != 0) {
                bail!("matrix.resolutions must start with the smallest positive resolution and be multiples of it")
            }
            let mut sorted = r.clone();
            sorted.sort_unstable();
            sorted.dedup();
            if &sorted != r {
                bail!("matrix.resolutions must be unique and ascending")
            }
        }
        if self.sam.min_mapq == 0 {
            bail!("sam.min_mapq must be positive")
        }
        for (label, p) in [
            ("read1", &self.read1),
            ("barcode_read", &self.barcode_read),
            ("read3", &self.read3),
            ("reference", &self.reference),
            ("chromsizes", &self.chromsizes),
            ("whitelist", &self.whitelist),
        ] {
            if !p.is_file() {
                bail!("{label} is not a regular file: {}", p.display())
            }
        }
        let suffixes: &[&str] = if self.aligner == "bwa" {
            &[".amb", ".ann", ".bwt", ".pac", ".sa"]
        } else {
            &[".amb", ".ann", ".bwt.2bit.64", ".pac"]
        };
        for suffix in suffixes {
            let mut p = self.aligner_index.clone().into_os_string();
            p.push(suffix);
            let p = PathBuf::from(p);
            if !p.is_file() {
                bail!("missing aligner index file: {}", p.display())
            }
        }
        if !self.aligner_executable.is_file() {
            bail!(
                "aligner executable is not a file: {}",
                self.aligner_executable.display()
            )
        }
        if !self.fastp_executable.is_file() {
            bail!(
                "fastp executable is not a file: {}",
                self.fastp_executable.display()
            )
        }
        Ok(())
    }

    pub(crate) fn scratch(&self) -> PathBuf {
        self.scratch_dir
            .clone()
            .unwrap_or_else(|| self.output_dir.join("scratch"))
    }
    pub(crate) fn index(&self) -> PathBuf {
        self.index_dir
            .clone()
            .unwrap_or_else(|| self.output_dir.join("fastq-index"))
    }
}
