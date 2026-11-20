//! Configuration for the droplet (3-read) Hi-C pipeline.
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashSet,
    fs::{self, File},
    path::{Path, PathBuf},
};

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct DropletConfig {
    /// Layout: `droplet` (3-read) or `two_read` (GAGE-seq / Spatial ATAC-Hi-C).
    #[serde(default = "default_layout")]
    pub sequencing_layout: String,
    /// R1: gDNA (Nextera read 1 side).
    pub read1: PathBuf,
    /// R2 for 3-read: cell barcode read.
    #[serde(default)]
    pub barcode_read: Option<PathBuf>,
    /// R3 for 3-read: gDNA mate.
    #[serde(default)]
    pub read3: Option<PathBuf>,
    /// R2 for 2-read: barcode rounds + linker + genomic mate.
    #[serde(default)]
    pub read2: Option<PathBuf>,
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
    /// 16 bp barcode whitelist (3-read; one barcode per line).
    #[serde(default)]
    pub whitelist: Option<PathBuf>,
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
    /// Target compressed chunk size accepted by the FASTQ indexer.
    #[serde(default = "default_chunk_size")]
    pub chunk_size: String,
    /// Boundary search window in MiB (increase for concatenated runs).
    #[serde(default = "default_search_mib")]
    pub index_search_mib: usize,
    /// Retain BAM (per-chunk BAMs are removed after a successful merge).
    #[serde(default = "no")]
    pub keep_bam: bool,
    /// Retain the final uncompressed pairs file (droplet default: true).
    #[serde(default = "yes")]
    pub keep_pairs: bool,
    /// Also write a BGZF-compressed `<prefix>.pairs.gz` with a pairix `.px2` index.
    #[serde(default)]
    pub compressed_pairs: bool,
    /// `fused` joins CB during parsing; `post` aligns first and joins afterwards.
    #[serde(default = "default_cb_join")]
    pub cb_join: String,
    /// In `post` mode, run the barcode prepass concurrently with alignment.
    #[serde(default)]
    pub prepass_overlap: bool,
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
    /// Two-read barcode rounds (takes precedence over `barcode` in two_read mode).
    #[serde(default)]
    pub cell_barcode: Option<CellBarcodeConfig>,
    /// Fixed 5'/3' trimming applied before alignment (two_read mode).
    #[serde(default)]
    pub additional_trimming: AdditionalTrimming,
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
    /// Decompression threads for the R2 barcode prepass (rapidgzip).
    #[serde(default = "default_barcode_threads")]
    pub threads: u32,
}

impl Default for BarcodeConfig {
    fn default() -> Self {
        Self {
            slice_start: 0,
            slice_len: 16,
            orientation: "rc".into(),
            max_mismatch: 1,
            threads: default_barcode_threads(),
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub(crate) struct AdditionalTrimming {
    pub read1_five_prime_bases: usize,
    pub read2_five_prime_bases: usize,
    pub read1_three_prime_bases: usize,
    pub read2_three_prime_bases: usize,
}

impl Default for AdditionalTrimming {
    fn default() -> Self {
        Self {
            read1_five_prime_bases: 0,
            read2_five_prime_bases: 0,
            read1_three_prime_bases: 0,
            read2_three_prime_bases: 0,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub(crate) struct CellBarcodeConfig {
    /// `forward`, `reverse_complement`/`rc`, or `auto`.
    pub barcode_orientation: String,
    pub barcode_rounds: Vec<BarcodeRound>,
    pub linker_regions: Vec<LinkerRegion>,
    pub barcode_structure: BarcodeStructure,
    pub valid_barcode_pairs_file: Option<PathBuf>,
    pub unique_molecule_identifier: Option<UniqueMoleculeIdentifier>,
}

impl Default for CellBarcodeConfig {
    fn default() -> Self {
        Self {
            barcode_orientation: "forward".into(),
            barcode_rounds: Vec::new(),
            linker_regions: Vec::new(),
            barcode_structure: BarcodeStructure::default(),
            valid_barcode_pairs_file: None,
            unique_molecule_identifier: None,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct BarcodeRound {
    pub round_name: String,
    #[serde(default)]
    pub first_base: Option<usize>,
    #[serde(default)]
    pub last_base: Option<usize>,
    #[serde(default)]
    pub segments: Vec<BarcodeSegment>,
    pub whitelist_file: PathBuf,
    #[serde(default)]
    pub whitelist_sequence_column: Option<usize>,
    #[serde(default)]
    pub whitelist_identifier_column: Option<usize>,
    #[serde(default = "default_identifier")]
    pub barcode_identifier: String,
    #[serde(default = "default_round_mismatch")]
    pub maximum_mismatches: u8,
    #[serde(default)]
    pub orientation: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct BarcodeSegment {
    pub first_base: usize,
    pub last_base: usize,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct LinkerRegion {
    pub region_name: String,
    pub first_base: usize,
    pub last_base: usize,
    #[serde(default)]
    pub expected_sequence: Option<String>,
    #[serde(default)]
    pub maximum_mismatches: usize,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub(crate) struct BarcodeStructure {
    pub order: Vec<String>,
    pub separator: String,
}

impl Default for BarcodeStructure {
    fn default() -> Self {
        Self {
            order: Vec::new(),
            separator: String::new(),
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub(crate) struct UniqueMoleculeIdentifier {
    pub enabled: bool,
    pub first_base: usize,
    pub last_base: usize,
}

impl Default for UniqueMoleculeIdentifier {
    fn default() -> Self {
        Self {
            enabled: false,
            first_base: 0,
            last_base: 0,
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

fn default_layout() -> String {
    "droplet".into()
}
fn default_identifier() -> String {
    "sequence".into()
}
fn default_round_mismatch() -> u8 {
    1
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
fn default_cb_join() -> String {
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
        c.read1 = resolve(&base, std::mem::take(&mut c.read1));
        c.reference = resolve(&base, std::mem::take(&mut c.reference));
        c.aligner_index = resolve(&base, std::mem::take(&mut c.aligner_index));
        c.aligner_executable = resolve(&base, std::mem::take(&mut c.aligner_executable));
        c.chromsizes = resolve(&base, std::mem::take(&mut c.chromsizes));
        c.output_dir = resolve(&base, std::mem::take(&mut c.output_dir));
        c.fastp_executable = resolve(&base, std::mem::take(&mut c.fastp_executable));
        for p in [
            &mut c.read2,
            &mut c.barcode_read,
            &mut c.read3,
            &mut c.whitelist,
        ] {
            if let Some(path) = p.as_mut() {
                *path = resolve(&base, std::mem::take(path));
            }
        }
        if let Some(cell) = c.cell_barcode.as_mut() {
            for round in &mut cell.barcode_rounds {
                round.whitelist_file =
                    resolve(&base, std::mem::take(&mut round.whitelist_file));
            }
            if let Some(path) = cell.valid_barcode_pairs_file.as_mut() {
                *path = resolve(&base, std::mem::take(path));
            }
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
        match self.sequencing_layout.as_str() {
            "droplet" | "two_read" => {}
            other => bail!("sequencing_layout must be droplet or two_read, got {other:?}"),
        }
        if self.is_two_read() {
            if self.barcode.slice_len == 0 || self.barcode.slice_len > 16 {
                // The legacy single-barcode settings are ignored in two_read mode.
            }
        } else {
            if self.barcode.slice_len == 0 || self.barcode.slice_len > 16 {
                bail!("barcode.slice_len must be 1..=16")
            }
            if self.barcode.max_mismatch > 2 {
                bail!("barcode.max_mismatch currently supports 0, 1 or 2")
            }
            if self.barcode.threads == 0 || self.barcode.threads > 32 {
                bail!("barcode.threads must be 1..=32")
            }
            match self.barcode.orientation.as_str() {
                "rc" | "forward" | "auto" => {}
                other => bail!(
                    "barcode.orientation must be rc, forward or auto, got {other:?}"
                ),
            }
        }
        match self.cb_join.as_str() {
            "fused" | "post" => {}
            other => bail!("cb_join must be fused or post, got {other:?}"),
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
            ("reference", &self.reference),
            ("chromsizes", &self.chromsizes),
        ] {
            if !p.is_file() {
                bail!("{label} is not a regular file: {}", p.display())
            }
        }
        if self.is_two_read() {
            let read2 = self
                .read2
                .as_deref()
                .context("read2 is required when sequencing_layout is two_read")?;
            if !read2.is_file() {
                bail!("read2 is not a regular file: {}", read2.display())
            }
            let cell = self
                .cell_barcode
                .as_ref()
                .context("cell_barcode is required when sequencing_layout is two_read")?;
            self.validate_two_read(cell)?;
        } else {
            for (label, p) in [
                ("barcode_read", &self.barcode_read),
                ("read3", &self.read3),
                ("whitelist", &self.whitelist),
            ] {
                let p = p
                    .as_deref()
                    .with_context(|| format!("{label} is required for droplet layout"))?;
                if !p.is_file() {
                    bail!("{label} is not a regular file: {}", p.display())
                }
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

    pub(crate) fn is_two_read(&self) -> bool {
        self.sequencing_layout == "two_read"
    }

    /// Second indexed mate: R3 for droplet, R2 for two_read.
    pub(crate) fn mate2(&self) -> Result<&Path> {
        if self.is_two_read() {
            self.read2
                .as_deref()
                .context("read2 is required when sequencing_layout is two_read")
        } else {
            self.read3
                .as_deref()
                .context("read3 is required for droplet layout")
        }
    }

    fn validate_two_read(&self, cell: &CellBarcodeConfig) -> Result<()> {
        if cell.barcode_rounds.is_empty() {
            bail!("cell_barcode.barcode_rounds must not be empty")
        }
        match cell.barcode_orientation.as_str() {
            "forward" | "rc" | "reverse_complement" | "auto" => {}
            other => bail!(
                "cell_barcode.barcode_orientation must be forward, rc or auto, got {other:?}"
            ),
        }
        let mut names = HashSet::new();
        for round in &cell.barcode_rounds {
            if round.round_name.is_empty() || !names.insert(round.round_name.clone()) {
                bail!("barcode round names must be non-empty and unique")
            }
            let segments: Vec<(usize, usize)> = if round.segments.is_empty() {
                vec![(
                    round
                        .first_base
                        .context("barcode round needs first_base/last_base or segments")?,
                    round
                        .last_base
                        .context("barcode round needs first_base/last_base or segments")?,
                )]
            } else {
                round
                    .segments
                    .iter()
                    .map(|s| (s.first_base, s.last_base))
                    .collect()
            };
            let mut previous_end = 0;
            for (first, last) in segments {
                if first == 0 || last < first || first <= previous_end {
                    bail!(
                        "round {} has invalid or unsorted 1-based segment {first}..{last}",
                        round.round_name
                    )
                }
                previous_end = last;
            }
            if !round.whitelist_file.is_file() {
                bail!(
                    "round {} whitelist is not a file: {}",
                    round.round_name,
                    round.whitelist_file.display()
                )
            }
            match round.barcode_identifier.as_str() {
                "sequence" | "well_name" => {}
                other => bail!(
                    "round {} barcode_identifier must be sequence or well_name, got {other:?}",
                    round.round_name
                ),
            }
            if round.maximum_mismatches > 2 {
                bail!(
                    "round {} maximum_mismatches must be 0..=2",
                    round.round_name
                )
            }
            if let Some(column) = round.whitelist_sequence_column {
                if column == 0 {
                    bail!("whitelist_sequence_column is 1-based")
                }
            }
            if let Some(column) = round.whitelist_identifier_column {
                if column == 0 {
                    bail!("whitelist_identifier_column is 1-based")
                }
            }
            if let Some(orientation) = &round.orientation {
                match orientation.as_str() {
                    "forward" | "rc" | "reverse_complement" => {}
                    other => bail!(
                        "round {} orientation must be forward or rc, got {other:?}",
                        round.round_name
                    ),
                }
            }
        }
        if !cell.barcode_structure.order.is_empty() {
            if cell.barcode_structure.order.len() != cell.barcode_rounds.len() {
                bail!("barcode_structure.order must list every barcode round")
            }
            for name in &cell.barcode_structure.order {
                if !names.contains(name) {
                    bail!("barcode_structure.order references unknown round {name:?}")
                }
            }
        }
        for linker in &cell.linker_regions {
            if linker.first_base == 0 || linker.last_base < linker.first_base {
                bail!(
                    "linker {} has invalid 1-based range",
                    linker.region_name
                )
            }
            let len = linker.last_base - linker.first_base + 1;
            if let Some(expected) = &linker.expected_sequence {
                if expected.len() != len {
                    bail!(
                        "linker {} expected_sequence length {} does not match region length {len}",
                        linker.region_name,
                        expected.len()
                    )
                }
                if !expected
                    .bytes()
                    .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'))
                {
                    bail!("linker {} expected_sequence must be ACGTN", linker.region_name)
                }
            }
            if linker.maximum_mismatches > len {
                bail!(
                    "linker {} maximum_mismatches exceeds region length",
                    linker.region_name
                )
            }
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
