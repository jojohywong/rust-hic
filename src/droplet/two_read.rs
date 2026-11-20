//! Two-read layout support (GAGE-seq / Spatial ATAC-Hi-C).
//!
//! R2 contains one or more barcode rounds, optional linker/detection regions,
//! and the genomic mate. The fused feeder matches the barcode rounds, stamps
//! `CB=<composite>;` onto both read names, and streams R1/R2 to fastp.
use super::config::{BarcodeRound, CellBarcodeConfig, DropletConfig};
use crate::fastq_index::PairedFastqReader;
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs,
    io::{BufWriter, Write},
    path::Path,
    sync::atomic::{AtomicBool, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};

const ORIENTATION_FILE: &str = "barcode_orientation_decision.yaml";
const ORIENTATION_SAMPLE: u64 = 50_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Status {
    Exact,
    Corrected,
    Ambiguous,
    Unassigned,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub(crate) struct RoundStats {
    pub exact: u64,
    pub corrected: u64,
    pub ambiguous: u64,
    pub unassigned: u64,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub(crate) struct LinkerStats {
    pub pass: u64,
    pub fail: u64,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub(crate) struct Stats {
    pub total: u64,
    pub assigned: u64,
    pub unassigned: u64,
    pub ambiguous: u64,
    pub short: u64,
    pub rounds: BTreeMap<String, RoundStats>,
    pub linkers: BTreeMap<String, LinkerStats>,
}

impl Stats {
    fn add(&mut self, other: &Stats) {
        self.total += other.total;
        self.assigned += other.assigned;
        self.unassigned += other.unassigned;
        self.ambiguous += other.ambiguous;
        self.short += other.short;
        for (name, s) in &other.rounds {
            let e = self.rounds.entry(name.clone()).or_default();
            e.exact += s.exact;
            e.corrected += s.corrected;
            e.ambiguous += s.ambiguous;
            e.unassigned += s.unassigned;
        }
        for (name, s) in &other.linkers {
            let e = self.linkers.entry(name.clone()).or_default();
            e.pass += s.pass;
            e.fail += s.fail;
        }
    }

    pub(crate) fn merge(&mut self, other: &Stats) {
        self.add(other)
    }

    pub(crate) fn assignable_fraction(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            self.assigned as f64 / self.total as f64
        }
    }
}

struct Entry {
    identifier: String,
    sequence: Vec<u8>,
}

struct Whitelist {
    entries: Vec<Entry>,
    by_sequence: HashMap<Vec<u8>, usize>,
    len: usize,
}

impl Whitelist {
    fn load(
        path: &Path,
        sequence_column: Option<usize>,
        identifier_column: Option<usize>,
    ) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("read whitelist {}", path.display()))?;
        let mut entries: Vec<Entry> = Vec::new();
        let mut by_sequence = HashMap::new();
        let mut len = None;
        for (line_no, line) in text.lines().enumerate() {
            let line = line.trim_end_matches('\r');
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (sequence, identifier) = if let Some(column) = sequence_column {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() < column {
                    bail!(
                        "{}:{}: whitelist has fewer than {} columns",
                        path.display(),
                        line_no + 1,
                        column
                    )
                }
                let sequence = fields[column - 1].trim();
                let identifier = identifier_column
                    .and_then(|c| fields.get(c - 1).map(|s| s.trim()))
                    .filter(|s| !s.is_empty())
                    .unwrap_or(sequence);
                (sequence.to_owned(), identifier.to_owned())
            } else {
                (line.trim().to_owned(), line.trim().to_owned())
            };
            if sequence.is_empty()
                || !sequence
                    .bytes()
                    .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
            {
                bail!("{}:{}: invalid barcode sequence", path.display(), line_no + 1)
            }
            match len {
                None => len = Some(sequence.len()),
                Some(n) if n != sequence.len() => bail!(
                    "{}:{}: barcode lengths differ ({} vs {})",
                    path.display(),
                    line_no + 1,
                    sequence.len(),
                    n
                ),
                _ => {}
            }
            let index = entries.len();
            let key = sequence.as_bytes().to_vec();
            if by_sequence.insert(key, index).is_some() {
                bail!("{}:{}: duplicate barcode sequence", path.display(), line_no + 1)
            }
            entries.push(Entry {
                identifier,
                sequence: sequence.into_bytes(),
            });
        }
        if entries.is_empty() {
            bail!("whitelist is empty: {}", path.display())
        }
        Ok(Self {
            entries,
            by_sequence,
            len: len.unwrap_or(0),
        })
    }

    fn classify(&self, sequence: &[u8], max_mismatches: u8) -> (Status, Option<usize>) {
        if sequence.len() != self.len {
            return (Status::Unassigned, None);
        }
        if let Some(&index) = self.by_sequence.get(sequence) {
            return (Status::Exact, Some(index));
        }
        let mut hit: Option<usize> = None;
        let mut buffer = sequence.to_vec();
        for position in 0..self.len {
            let original = buffer[position];
            for base in b"ACGT" {
                if *base == original {
                    continue;
                }
                buffer[position] = *base;
                if let Some(&index) = self.by_sequence.get(&buffer) {
                    if hit.is_some_and(|h| h != index) {
                        return (Status::Ambiguous, None);
                    }
                    hit = Some(index);
                }
            }
            buffer[position] = original;
        }
        if let Some(index) = hit {
            return (Status::Corrected, Some(index));
        }
        if max_mismatches >= 2 {
            for first in 0..self.len {
                let original_first = buffer[first];
                for base_first in b"ACGT" {
                    if *base_first == original_first {
                        continue;
                    }
                    buffer[first] = *base_first;
                    for second in first + 1..self.len {
                        let original_second = buffer[second];
                        for base_second in b"ACGT" {
                            if *base_second == original_second {
                                continue;
                            }
                            buffer[second] = *base_second;
                            if let Some(&index) = self.by_sequence.get(&buffer) {
                                if hit.is_some_and(|h| h != index) {
                                    buffer[second] = original_second;
                                    buffer[first] = original_first;
                                    return (Status::Ambiguous, None);
                                }
                                hit = Some(index);
                            }
                        }
                        buffer[second] = original_second;
                    }
                    buffer[first] = original_first;
                }
            }
        }
        match hit {
            Some(index) => (Status::Corrected, Some(index)),
            None => (Status::Unassigned, None),
        }
    }
}

struct Round {
    name: String,
    segments: Vec<(usize, usize)>, // 0-based half-open
    whitelist: Whitelist,
    identifier_kind: IdentifierKind,
    max_mismatches: u8,
    reverse_complement: bool,
    end: usize,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum IdentifierKind {
    Sequence,
    WellName,
}

struct Linker {
    name: String,
    start: usize,
    end: usize,
    expected: Option<Vec<u8>>,
    max_mismatches: usize,
}

pub(crate) struct Call {
    pub(crate) cb: Option<String>,
    pub(crate) rounds: Vec<(String, Status)>,
    pub(crate) linkers: Vec<(String, bool)>,
    pub(crate) short: bool,
    pub(crate) exact_all: bool,
}

pub(crate) struct Matcher {
    rounds: Vec<Round>,
    order: Vec<usize>,
    separator: String,
    linkers: Vec<Linker>,
    valid_pairs: Option<HashSet<String>>,
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence.iter().rev().map(|b| complement(*b)).collect()
}

fn normalize_orientation(value: &str) -> Result<bool> {
    match value {
        "forward" => Ok(false),
        "rc" | "reverse_complement" => Ok(true),
        other => bail!("barcode orientation must be forward or rc, got {other:?}"),
    }
}

impl Matcher {
    pub(crate) fn load(c: &DropletConfig, orientation: bool) -> Result<Self> {
        let cell = c
            .cell_barcode
            .as_ref()
            .context("cell_barcode is required for two_read")?;
        let mut rounds = Vec::new();
        let mut by_name = HashMap::new();
        for round in &cell.barcode_rounds {
            let (segments, end) = round_segments(round)?;
            let whitelist = Whitelist::load(
                &round.whitelist_file,
                round.whitelist_sequence_column,
                round.whitelist_identifier_column,
            )?;
            let identifier_kind = match round.barcode_identifier.as_str() {
                "sequence" => IdentifierKind::Sequence,
                "well_name" => IdentifierKind::WellName,
                other => bail!("barcode_identifier must be sequence or well_name, got {other:?}"),
            };
            let reverse_complement = match round.orientation.as_deref() {
                Some(value) => normalize_orientation(value)?,
                None => orientation,
            };
            by_name.insert(round.round_name.clone(), rounds.len());
            rounds.push(Round {
                name: round.round_name.clone(),
                segments,
                whitelist,
                identifier_kind,
                max_mismatches: round.maximum_mismatches,
                reverse_complement,
                end,
            });
        }
        let order = if cell.barcode_structure.order.is_empty() {
            (0..rounds.len()).collect()
        } else {
            cell.barcode_structure
                .order
                .iter()
                .map(|name| {
                    by_name
                        .get(name)
                        .copied()
                        .with_context(|| format!("unknown barcode round {name:?}"))
                })
                .collect::<Result<Vec<_>>>()?
        };
        let linkers = cell
            .linker_regions
            .iter()
            .map(|linker| {
                let start = linker.first_base - 1;
                let end = linker.last_base;
                let expected = linker
                    .expected_sequence
                    .as_ref()
                    .map(|s| s.as_bytes().to_vec());
                Ok(Linker {
                    name: linker.region_name.clone(),
                    start,
                    end,
                    expected,
                    max_mismatches: linker.maximum_mismatches,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let valid_pairs = cell
            .valid_barcode_pairs_file
            .as_ref()
            .map(|path| -> Result<HashSet<String>> {
                let text = fs::read_to_string(path)
                    .with_context(|| format!("read {}", path.display()))?;
                Ok(text
                    .lines()
                    .map(|line| line.trim().to_owned())
                    .filter(|line| !line.is_empty())
                    .collect())
            })
            .transpose()?;
        Ok(Self {
            rounds,
            order,
            separator: cell.barcode_structure.separator.clone(),
            linkers,
            valid_pairs,
        })
    }

    pub(crate) fn classify(&self, read2: &[u8]) -> Call {
        let mut short = false;
        let mut matched = vec![None; self.rounds.len()];
        let mut round_status = Vec::with_capacity(self.rounds.len());
        for (index, round) in self.rounds.iter().enumerate() {
            if round.end > read2.len() {
                short = true;
                round_status.push((round.name.clone(), Status::Unassigned));
                continue;
            }
            let mut sequence = Vec::new();
            for (start, end) in &round.segments {
                let piece = &read2[*start..*end];
                if round.reverse_complement {
                    sequence.extend(reverse_complement(piece));
                } else {
                    sequence.extend_from_slice(piece);
                }
            }
            let (status, hit) = round.whitelist.classify(&sequence, round.max_mismatches);
            matched[index] = hit;
            round_status.push((round.name.clone(), status));
        }
        let mut linkers = Vec::with_capacity(self.linkers.len());
        let mut linkers_ok = true;
        for linker in &self.linkers {
            let ok = if linker.end > read2.len() {
                short = true;
                false
            } else if let Some(expected) = &linker.expected {
                let observed = &read2[linker.start..linker.end];
                hamming_with_wildcards(expected, observed) <= linker.max_mismatches
            } else {
                true
            };
            linkers_ok &= ok;
            linkers.push((linker.name.clone(), ok));
        }
        let mut exact_all = linkers_ok && !short;
        for (round, status) in self.rounds.iter().zip(round_status.iter()) {
            if status.1 != Status::Exact {
                exact_all = false;
            }
            let _ = round;
        }
        let mut cb = None;
        if linkers_ok
            && round_status
                .iter()
                .all(|(_, status)| matches!(status, Status::Exact | Status::Corrected))
        {
            let mut parts = Vec::with_capacity(self.order.len());
            for &index in &self.order {
                let round = &self.rounds[index];
                let entry = &round.whitelist.entries[matched[index].expect("round matched")];
                parts.push(match round.identifier_kind {
                    IdentifierKind::Sequence => {
                        String::from_utf8_lossy(&entry.sequence).into_owned()
                    }
                    IdentifierKind::WellName => entry.identifier.clone(),
                });
            }
            let composite = parts.join(&self.separator);
            if self
                .valid_pairs
                .as_ref()
                .is_none_or(|pairs| pairs.contains(&composite))
            {
                cb = Some(composite);
            }
        }
        Call {
            cb,
            rounds: round_status,
            linkers,
            short,
            exact_all,
        }
    }
}

fn round_segments(round: &BarcodeRound) -> Result<(Vec<(usize, usize)>, usize)> {
    let config_segments: Vec<(usize, usize)> = if round.segments.is_empty() {
        vec![(
            round
                .first_base
                .context("round needs first_base/last_base or segments")?,
            round
                .last_base
                .context("round needs first_base/last_base or segments")?,
        )]
    } else {
        round
            .segments
            .iter()
            .map(|s| (s.first_base, s.last_base))
            .collect()
    };
    let mut segments = Vec::new();
    let mut end = 0;
    for (first, last) in config_segments {
        let start = first
            .checked_sub(1)
            .context("barcode positions are 1-based")?;
        if last < first {
            bail!("invalid barcode segment {first}..{last}")
        }
        end = end.max(last);
        segments.push((start, last));
    }
    Ok((segments, end))
}

fn hamming_with_wildcards(expected: &[u8], observed: &[u8]) -> usize {
    expected
        .iter()
        .zip(observed)
        .filter(|(e, o)| **e != b'N' && *e != *o)
        .count()
}

fn record_stats(stats: &mut Stats, call: &Call) {
    stats.total += 1;
    if call.short {
        stats.short += 1;
    }
    for (name, status) in &call.rounds {
        let round = stats.rounds.entry(name.clone()).or_default();
        match status {
            Status::Exact => round.exact += 1,
            Status::Corrected => round.corrected += 1,
            Status::Ambiguous => round.ambiguous += 1,
            Status::Unassigned => round.unassigned += 1,
        }
    }
    for (name, ok) in &call.linkers {
        let linker = stats.linkers.entry(name.clone()).or_default();
        if *ok {
            linker.pass += 1;
        } else {
            linker.fail += 1;
        }
    }
    if call.cb.is_some() {
        stats.assigned += 1;
    } else if call
        .rounds
        .iter()
        .any(|(_, status)| *status == Status::Ambiguous)
    {
        stats.ambiguous += 1;
    } else {
        stats.unassigned += 1;
    }
}

fn write_record<W: Write>(
    out: &mut W,
    cb: &str,
    name: &[u8],
    sequence: &[u8],
    quality: &[u8],
) -> Result<()> {
    out.write_all(b"@CB=")?;
    out.write_all(cb.as_bytes())?;
    out.write_all(b";")?;
    out.write_all(name)?;
    out.write_all(b"\n")?;
    out.write_all(sequence)?;
    out.write_all(b"\n+\n")?;
    out.write_all(quality)?;
    out.write_all(b"\n")?;
    Ok(())
}

pub(crate) fn feeder<W: Write>(
    c: &DropletConfig,
    index_dir: &Path,
    chunk: usize,
    output: W,
    cancel: &AtomicBool,
) -> Result<Stats> {
    let orientation = load_orientation(c)?;
    let matcher = Matcher::load(c, orientation)?;
    let mut reader = PairedFastqReader::open(index_dir, chunk)?;
    let mut out = BufWriter::with_capacity(1024 * 1024, output);
    let mut stats = Stats::default();
    while let Some(pair) = reader.next_pair()? {
        if cancel.load(Ordering::Relaxed) {
            bail!("FASTQ stream cancelled")
        }
        let call = matcher.classify(&pair.seq2);
        record_stats(&mut stats, &call);
        let cb = call.cb.as_deref().unwrap_or(".");
        write_record(&mut out, cb, &pair.name, &pair.seq1, &pair.qual1)?;
        write_record(&mut out, cb, &pair.name, &pair.seq2, &pair.qual2)?;
    }
    out.flush()?;
    Ok(stats)
}

#[derive(Serialize)]
struct DecisionFile {
    barcode_orientation_decision: Decision,
}

#[derive(Serialize, Deserialize)]
struct Decision {
    schema_version: u32,
    sequencing_layout: String,
    chosen_orientation: String,
    decision_rule: String,
    reason: String,
    sampled_reads: u64,
    candidates: BTreeMap<String, Candidate>,
    input_fingerprint: Fingerprint,
    decided_at_unix_seconds: u64,
}

#[derive(Serialize, Deserialize, Default)]
struct Candidate {
    reads_with_all_rounds_exact: u64,
    fraction: f64,
}

#[derive(Serialize, Deserialize, Default)]
struct Fingerprint {
    read2_fastq: String,
    read2_size_bytes: u64,
    read2_modified_unix_seconds: u64,
    whitelists: Vec<String>,
}

fn metadata_fingerprint(c: &DropletConfig, cell: &CellBarcodeConfig) -> Fingerprint {
    let read2 = c.read2.as_deref();
    let (path, size, modified) = read2
        .and_then(|path| fs::metadata(path).ok().map(|m| (path, m)))
        .map(|(path, m)| {
            (
                path.display().to_string(),
                m.len(),
                m.modified()
                    .ok()
                    .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
                    .map(|d| d.as_secs())
                    .unwrap_or(0),
            )
        })
        .unwrap_or_default();
    Fingerprint {
        read2_fastq: path,
        read2_size_bytes: size,
        read2_modified_unix_seconds: modified,
        whitelists: cell
            .barcode_rounds
            .iter()
            .map(|r| r.whitelist_file.display().to_string())
            .collect(),
    }
}

fn write_decision(
    c: &DropletConfig,
    cell: &CellBarcodeConfig,
    chosen: &str,
    rule: &str,
    reason: &str,
    sampled: u64,
    candidates: BTreeMap<String, Candidate>,
) -> Result<()> {
    let file = DecisionFile {
        barcode_orientation_decision: Decision {
            schema_version: 1,
            sequencing_layout: c.sequencing_layout.clone(),
            chosen_orientation: chosen.to_owned(),
            decision_rule: rule.to_owned(),
            reason: reason.to_owned(),
            sampled_reads: sampled,
            candidates,
            input_fingerprint: metadata_fingerprint(c, cell),
            decided_at_unix_seconds: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map(|d| d.as_secs())
                .unwrap_or(0),
        },
    };
    let path = c.output_dir.join(ORIENTATION_FILE);
    fs::write(&path, serde_yaml::to_string(&file)?)
        .with_context(|| format!("write {}", path.display()))
}

pub(crate) fn resolve_orientation(c: &DropletConfig, index_dir: &Path) -> Result<String> {
    let cell = c
        .cell_barcode
        .as_ref()
        .context("cell_barcode is required for two_read")?;
    if cell.barcode_orientation != "auto" {
        let rc = normalize_orientation(&cell.barcode_orientation)?;
        let chosen = if rc { "reverse_complement" } else { "forward" };
        write_decision(
            c,
            cell,
            chosen,
            "explicit in config",
            "barcode_orientation is explicit",
            0,
            BTreeMap::new(),
        )?;
        return Ok(chosen.to_owned());
    }
    let forward = Matcher::load(c, false)?;
    let reverse = Matcher::load(c, true)?;
    let mut reader = PairedFastqReader::open(index_dir, 0)?;
    let (mut fw, mut rc, mut seen) = (0u64, 0u64, 0u64);
    while seen < ORIENTATION_SAMPLE {
        let Some(pair) = reader.next_pair()? else {
            break;
        };
        seen += 1;
        if forward.classify(&pair.seq2).exact_all {
            fw += 1;
        }
        if reverse.classify(&pair.seq2).exact_all {
            rc += 1;
        }
    }
    if seen == 0 {
        bail!("cannot detect barcode orientation: R2 is empty")
    }
    let mut candidates = BTreeMap::new();
    candidates.insert(
        "forward".to_owned(),
        Candidate {
            reads_with_all_rounds_exact: fw,
            fraction: fw as f64 / seen as f64,
        },
    );
    candidates.insert(
        "reverse_complement".to_owned(),
        Candidate {
            reads_with_all_rounds_exact: rc,
            fraction: rc as f64 / seen as f64,
        },
    );
    let (best, best_count, other) = if fw >= rc {
        ("forward", fw, rc)
    } else {
        ("reverse_complement", rc, fw)
    };
    let clear = best_count >= 5 * other.max(1) && best_count * 2 >= seen;
    if !clear {
        write_decision(
            c,
            cell,
            "undetermined",
            "best >= 5x runner-up and >= 50% of sampled reads",
            "orientation is ambiguous; set barcode_orientation explicitly",
            seen,
            candidates,
        )?;
        bail!(
            "could not determine barcode orientation clearly (forward {fw}, reverse {rc} of {seen}); set barcode_orientation explicitly"
        )
    }
    write_decision(
        c,
        cell,
        best,
        "best >= 5x runner-up and >= 50% of sampled reads",
        &format!("{best} won clearly"),
        seen,
        candidates,
    )?;
    Ok(best.to_owned())
}

pub(crate) fn load_orientation(c: &DropletConfig) -> Result<bool> {
    let cell = c
        .cell_barcode
        .as_ref()
        .context("cell_barcode is required for two_read")?;
    if cell.barcode_orientation != "auto" {
        return normalize_orientation(&cell.barcode_orientation);
    }
    let path = c.output_dir.join(ORIENTATION_FILE);
    let text = fs::read_to_string(&path).with_context(|| {
        format!(
            "barcode orientation is auto but {} is missing; run the prepare/orientation step first",
            path.display()
        )
    })?;
    #[derive(Deserialize)]
    struct File {
        barcode_orientation_decision: Body,
    }
    #[derive(Deserialize)]
    struct Body {
        chosen_orientation: String,
    }
    let file: File = serde_yaml::from_str(&text)
        .with_context(|| format!("parse {}", path.display()))?;
    normalize_orientation(&file.barcode_orientation_decision.chosen_orientation)
}

pub(crate) fn write_qc(c: &DropletConfig, qc_dir: &Path) -> Result<()> {
    let mut total = Stats::default();
    let mut per_chunk = Vec::new();
    let mut entries: Vec<_> = fs::read_dir(qc_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("barcode-chunk-") && n.ends_with(".json"))
        })
        .collect();
    entries.sort();
    for path in entries {
        let stats: Stats = serde_json::from_reader(
            fs::File::open(&path).with_context(|| format!("open {}", path.display()))?,
        )
        .with_context(|| format!("parse {}", path.display()))?;
        total.merge(&stats);
        let chunk = path
            .file_stem()
            .and_then(|s| s.to_str())
            .and_then(|s| s.strip_prefix("barcode-chunk-"))
            .unwrap_or("?")
            .to_owned();
        per_chunk.push((chunk, stats));
    }
    let mut out = BufWriter::new(fs::File::create(qc_dir.join("barcode_stats.tsv"))?);
    let orientation = load_orientation(c).unwrap_or(false);
    writeln!(
        out,
        "#orientation: {}",
        if orientation {
            "reverse_complement"
        } else {
            "forward"
        }
    )?;
    writeln!(
        out,
        "chunk\ttotal\tassigned\tunassigned\tambiguous\tshort\tassignable_fraction"
    )?;
    for (chunk, stats) in &per_chunk {
        writeln!(
            out,
            "{chunk}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
            stats.total,
            stats.assigned,
            stats.unassigned,
            stats.ambiguous,
            stats.short,
            stats.assignable_fraction()
        )?;
    }
    writeln!(
        out,
        "total\t{}\t{}\t{}\t{}\t{}\t{:.6}",
        total.total,
        total.assigned,
        total.unassigned,
        total.ambiguous,
        total.short,
        total.assignable_fraction()
    )?;
    let mut rounds = BufWriter::new(fs::File::create(qc_dir.join("barcode_round_stats.tsv"))?);
    writeln!(rounds, "chunk\tround\ttotal\texact\tcorrected\tambiguous\tunassigned")?;
    for (chunk, stats) in &per_chunk {
        for (round, s) in &stats.rounds {
            let n = s.exact + s.corrected + s.ambiguous + s.unassigned;
            writeln!(
                rounds,
                "{chunk}\t{round}\t{n}\t{}\t{}\t{}\t{}",
                s.exact, s.corrected, s.ambiguous, s.unassigned
            )?;
        }
    }
    let mut linkers =
        BufWriter::new(fs::File::create(qc_dir.join("linker_region_stats.tsv"))?);
    writeln!(linkers, "chunk\tregion\tpass\tfail")?;
    for (chunk, stats) in &per_chunk {
        for (linker, s) in &stats.linkers {
            writeln!(linkers, "{chunk}\t{linker}\t{}\t{}", s.pass, s.fail)?;
        }
    }
    let _ = c;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wildcard_linker_matching() {
        assert_eq!(hamming_with_wildcards(b"TTNA", b"TTCA"), 0);
        assert_eq!(hamming_with_wildcards(b"TTNA", b"TTGA"), 0);
        assert_eq!(hamming_with_wildcards(b"TTNA", b"TAGA"), 1);
    }

    #[test]
    fn rc_helper() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT".to_vec());
        assert_eq!(reverse_complement(b"AACC"), b"GGTT".to_vec());
    }
}
