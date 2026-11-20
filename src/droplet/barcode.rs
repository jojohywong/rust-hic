//! Dedicated R2 barcode worker: one sequential pass over the barcode read.
//!
//! R2 is not part of the R1/R3 FASTQ index. Since all three files are
//! record-synchronized, the R1/R3 boundary read IDs delimit R2 chunks too.
//! The pass writes one small `readID<TAB>CB` sidecar per chunk; chunk workers
//! join CB to SAM records by QNAME.
use crate::fastq_index::rapidgzip::RapidgzipReader;
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    cell::{Cell, RefCell},
    collections::{HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use super::{config::DropletConfig, index::IndexManifest};

/// Bump when barcode matching changes; prevents reuse of stale sidecars.
const MATCHER_VERSION: u32 = 2;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub(crate) struct BarcodeStats {
    pub total: u64,
    pub exact: u64,
    pub corrected: u64,
    pub ambiguous: u64,
    pub poly_g: u64,
    pub unassigned: u64,
}

impl BarcodeStats {
    fn add(&mut self, other: &Self) {
        self.total += other.total;
        self.exact += other.exact;
        self.corrected += other.corrected;
        self.ambiguous += other.ambiguous;
        self.poly_g += other.poly_g;
        self.unassigned += other.unassigned;
    }
}

pub(crate) struct BarcodePreprocessing {
    #[allow(dead_code)]
    pub dir: PathBuf,
    #[allow(dead_code)]
    pub chunks: usize,
    pub stats: Vec<BarcodeStats>,
    pub total: BarcodeStats,
    pub orientation: &'static str,
}

impl BarcodePreprocessing {
    #[allow(dead_code)]
    pub(crate) fn sidecar(&self, chunk: usize) -> PathBuf {
        self.dir.join(format!("chunk-{chunk:04}.tsv"))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Match {
    Exact([u8; 16]),
    Corrected([u8; 16]),
    Ambiguous,
    PolyG,
    Unassigned,
}

struct Matcher {
    whitelist: HashSet<[u8; 16]>,
    entries: Vec<[u8; 16]>,
    /// Index of the three barcode thirds, used for a pigeonhole search of
    /// neighbors within two mismatches (at least one third must match exactly).
    part_a: HashMap<[u8; 6], Vec<u32>>,
    part_b: HashMap<[u8; 6], Vec<u32>>,
    part_c: HashMap<[u8; 6], Vec<u32>>,
    bounds: (usize, usize),
    start: usize,
    len: usize,
    max_mismatch: u8,
    rc: bool,
}

thread_local! {
    /// Per-thread dedup marks for the two-mismatch candidate search.
    static SEEN: RefCell<Vec<u32>> = const { RefCell::new(Vec::new()) };
    static GENERATION: Cell<u32> = const { Cell::new(0) };
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

fn reverse_complement(seq: &[u8], out: &mut [u8; 16]) {
    for (i, b) in seq.iter().rev().enumerate() {
        out[i] = complement(*b);
    }
}

fn part_key(seq: &[u8]) -> [u8; 6] {
    let mut k = [0u8; 6];
    k[..seq.len()].copy_from_slice(seq);
    k
}

impl Matcher {
    fn load(c: &DropletConfig, orientation: &'static str) -> Result<Self> {
        let whitelist = c
            .whitelist
            .as_deref()
            .context("whitelist is required for the droplet barcode prepass")?;
        let file = File::open(whitelist)
            .with_context(|| format!("open whitelist {}", whitelist.display()))?;
        let mut entries = Vec::new();
        for (i, line) in BufReader::new(file).lines().enumerate() {
            let line = line?;
            let b = line.trim();
            if b.is_empty() || b.starts_with('#') {
                continue;
            }
            if b.len() != c.barcode.slice_len {
                bail!(
                    "whitelist line {} has length {} but barcode.slice_len is {}",
                    i + 1,
                    b.len(),
                    c.barcode.slice_len
                )
            }
            if b.len() > 16 {
                bail!("barcode.slice_len > 16 is not supported")
            }
            if !b.bytes().all(|x| matches!(x, b'A' | b'C' | b'G' | b'T')) {
                bail!("whitelist line {} contains non-ACGT bases", i + 1)
            }
            let mut a = [0u8; 16];
            a[..b.len()].copy_from_slice(b.as_bytes());
            entries.push(a);
        }
        if entries.is_empty() {
            bail!("whitelist is empty: {}", whitelist.display())
        }
        Ok(Self::from_entries(
            entries,
            c.barcode.slice_start,
            c.barcode.slice_len,
            c.barcode.max_mismatch,
            orientation == "rc",
        ))
    }

    fn from_entries(
        entries: Vec<[u8; 16]>,
        start: usize,
        len: usize,
        max_mismatch: u8,
        rc: bool,
    ) -> Self {
        let whitelist: HashSet<[u8; 16]> = entries.iter().copied().collect();
        // Three parts with two mismatches: pigeonhole guarantees that at least
        // one part is an exact match, so a part lookup cannot miss a neighbor.
        let a = len / 3;
        let b = (2 * len) / 3;
        let mut part_a: HashMap<[u8; 6], Vec<u32>> = HashMap::new();
        let mut part_b: HashMap<[u8; 6], Vec<u32>> = HashMap::new();
        let mut part_c: HashMap<[u8; 6], Vec<u32>> = HashMap::new();
        for (i, e) in entries.iter().enumerate() {
            part_a.entry(part_key(&e[..a])).or_default().push(i as u32);
            part_b.entry(part_key(&e[a..b])).or_default().push(i as u32);
            part_c.entry(part_key(&e[b..len])).or_default().push(i as u32);
        }
        Self {
            whitelist,
            entries,
            part_a,
            part_b,
            part_c,
            bounds: (a, b),
            start,
            len,
            max_mismatch,
            rc,
        }
    }

    fn classify(&self, seq: &[u8]) -> Match {
        let Some(slice) = seq.get(self.start..self.start + self.len) else {
            return Match::Unassigned;
        };
        let mut bc = [0u8; 16];
        if self.rc {
            reverse_complement(slice, &mut bc);
        } else {
            bc[..self.len].copy_from_slice(slice);
        }
        // A read made only of G (and N) is the classic no-barcode poly-G
        // artifact. Individual N bases are *not* poly-G: they are treated as
        // unknown bases and may still be corrected against the whitelist.
        if bc[..self.len]
            .iter()
            .all(|x| matches!(x, b'G' | b'N'))
        {
            return Match::PolyG;
        }
        let non_acgt = bc[..self.len]
            .iter()
            .filter(|x| !matches!(x, b'A' | b'C' | b'G' | b'T'))
            .count();
        if non_acgt > self.max_mismatch as usize {
            return Match::Unassigned;
        }
        let mut key = [0u8; 16];
        key[..self.len].copy_from_slice(&bc[..self.len]);
        if self.whitelist.contains(&key) {
            return Match::Exact(key);
        }
        if self.max_mismatch == 0 {
            return Match::Unassigned;
        }
        // One substitution: this also resolves single N bases, since N is
        // replaced by each of the four bases in turn.
        let mut hits = 0usize;
        let mut first = key;
        'outer: for pos in 0..self.len {
            let original = key[pos];
            for alt in b"ACGT" {
                if *alt == original {
                    continue;
                }
                key[pos] = *alt;
                if self.whitelist.contains(&key) {
                    hits += 1;
                    if hits == 1 {
                        first = key;
                    } else {
                        break 'outer;
                    }
                }
            }
            key[pos] = original;
        }
        match hits {
            0 => {}
            1 => return Match::Corrected(first),
            _ => return Match::Ambiguous,
        }
        if self.max_mismatch < 2 {
            return Match::Unassigned;
        }
        // Two substitutions via the three-part index.
        let (a, b) = self.bounds;
        let ka = part_key(&bc[..a]);
        let kb = part_key(&bc[a..b]);
        let kc = part_key(&bc[b..self.len]);
        let generation = GENERATION.with(|g| {
            let next = g.get().wrapping_add(1);
            let next = if next == 0 { 1 } else { next };
            g.set(next);
            next
        });
        let (hits, first) = SEEN.with(|cell| {
            let mut seen = cell.borrow_mut();
            if seen.len() < self.entries.len() {
                seen.resize(self.entries.len(), 0);
            }
            if generation == 1 {
                seen.fill(0);
            }
            let mut hits = 0usize;
            let mut first = [0u8; 16];
            for list in [
                self.part_a.get(&ka),
                self.part_b.get(&kb),
                self.part_c.get(&kc),
            ]
            .into_iter()
            .flatten()
            {
                for &i in list {
                    let idx = i as usize;
                    if seen[idx] == generation {
                        continue;
                    }
                    seen[idx] = generation;
                    let cand = &self.entries[idx];
                    let mut dist = 0usize;
                    for pos in 0..self.len {
                        if bc[pos] != cand[pos] {
                            dist += 1;
                            if dist > 2 {
                                break;
                            }
                        }
                    }
                    if dist <= 2 {
                        hits += 1;
                        if hits == 1 {
                            first = *cand;
                        } else {
                            return (2usize, first);
                        }
                    }
                }
            }
            (hits, first)
        });
        if hits == 1 {
            Match::Corrected(first)
        } else if hits > 1 {
            Match::Ambiguous
        } else {
            Match::Unassigned
        }
    }
}

fn trim_newline(line: &mut Vec<u8>) {
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
}

fn name_of(header: &[u8]) -> &[u8] {
    let end = header
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(header.len());
    &header[..end]
}

/// Open R2 with rapidgzip (parallel) for gzip input, or plain buffered I/O.
fn open_barcode_read(path: &Path, threads: u32) -> Result<Box<dyn BufRead>> {
    let capacity = 1024 * 1024;
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        let reader = RapidgzipReader::open(path, threads)
            .with_context(|| format!("open barcode read {}", path.display()))?;
        Ok(Box::new(BufReader::with_capacity(capacity, reader)))
    } else {
        let file =
            File::open(path).with_context(|| format!("open barcode read {}", path.display()))?;
        Ok(Box::new(BufReader::with_capacity(capacity, file)))
    }
}

/// Guess orientation on a small prefix of R2 by counting exact matches.
fn detect_orientation(c: &DropletConfig) -> Result<&'static str> {
    let barcode_read = c
        .barcode_read
        .as_deref()
        .context("barcode_read is required for the droplet barcode preprocessing")?;
    let mut reader = open_barcode_read(barcode_read, c.barcode.matcher_threads)?;
    let (matcher_rc, matcher_fw) = (
        Matcher::load(c, "rc")?,
        Matcher::load(c, "forward")?,
    );
    let (mut rc, mut fw, mut seen) = (0u64, 0u64, 0u64);
    let (mut h, mut s, mut p, mut q) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    while seen < 20_000 {
        h.clear();
        s.clear();
        p.clear();
        q.clear();
        if reader.read_until(b'\n', &mut h)? == 0
            || reader.read_until(b'\n', &mut s)? == 0
            || reader.read_until(b'\n', &mut p)? == 0
            || reader.read_until(b'\n', &mut q)? == 0
        {
            break;
        }
        trim_newline(&mut s);
        if matches!(
            matcher_rc.classify(&s),
            Match::Exact(_) | Match::Corrected(_)
        ) {
            rc += 1;
        }
        if matches!(
            matcher_fw.classify(&s),
            Match::Exact(_) | Match::Corrected(_)
        ) {
            fw += 1;
        }
        seen += 1;
    }
    if rc == 0 && fw == 0 {
        bail!("could not detect barcode orientation: no exact matches in first 20,000 R2 reads")
    }
    Ok(if rc >= fw { "rc" } else { "forward" })
}

const PREPASS_BLOCK: usize = 65_536;

#[allow(clippy::too_many_arguments)]
fn flush_block(
    pool: &rayon::ThreadPool,
    matcher: &Matcher,
    data: &mut Vec<u8>,
    spans: &mut Vec<(u32, u32, u32, u32)>,
    stats: &mut BarcodeStats,
    writer: &mut BufWriter<File>,
) -> Result<usize> {
    if spans.is_empty() {
        return Ok(0);
    }
    let matches: Vec<Match> = pool.install(|| {
        spans
            .par_iter()
            .map(|&(_, _, ss, se)| matcher.classify(&data[ss as usize..se as usize]))
            .collect()
    });
    for (i, &(ns, ne, _, _)) in spans.iter().enumerate() {
        let m = matches[i];
        stats.total += 1;
        let mut cb_buf = [0u8; 16];
        let cb_len = match m {
            Match::Exact(b) | Match::Corrected(b) => {
                if matches!(m, Match::Exact(_)) {
                    stats.exact += 1;
                } else {
                    stats.corrected += 1;
                }
                cb_buf = b;
                matcher.len
            }
            Match::Ambiguous => {
                stats.ambiguous += 1;
                cb_buf[0] = b'.';
                1
            }
            Match::PolyG => {
                stats.poly_g += 1;
                cb_buf[0] = b'.';
                1
            }
            Match::Unassigned => {
                stats.unassigned += 1;
                cb_buf[0] = b'.';
                1
            }
        };
        writer.write_all(&data[ns as usize..ne as usize])?;
        writer.write_all(b"\t")?;
        writer.write_all(&cb_buf[..cb_len])?;
        writer.write_all(b"\n")?;
    }
    let n = spans.len();
    data.clear();
    spans.clear();
    Ok(n)
}

pub(crate) fn preprocess(c: &DropletConfig, index: &IndexManifest, dir: &Path) -> Result<BarcodePreprocessing> {
    let stats_path = dir.join("stats.json");
    if c.reuse_barcodes && stats_path.is_file() {
        #[derive(Deserialize)]
        struct Saved {
            orientation: String,
            stats: Vec<BarcodeStats>,
            #[serde(default)]
            slice_start: usize,
            #[serde(default)]
            slice_len: usize,
            #[serde(default)]
            max_mismatch: u8,
            #[serde(default)]
            matcher_version: u32,
        }
        let saved: Saved = serde_json::from_reader(File::open(&stats_path)?)?;
        let compatible = saved.slice_start == c.barcode.slice_start
            && saved.slice_len == c.barcode.slice_len
            && saved.max_mismatch == c.barcode.max_mismatch
            && saved.matcher_version == MATCHER_VERSION;
        if compatible && saved.stats.len() == index.chunks {
            for i in 0..index.chunks {
                if !dir.join(format!("chunk-{i:04}.tsv")).is_file() {
                    bail!("barcode sidecar missing for chunk {i}; rerun without reuse_barcodes")
                }
            }
            let orientation = match saved.orientation.as_str() {
                "rc" => "rc",
                "forward" => "forward",
                other => bail!("invalid saved orientation {other:?}"),
            };
            let mut total = BarcodeStats::default();
            for s in &saved.stats {
                total.add(s);
            }
            eprintln!("event=barcode_reuse chunks={}", index.chunks);
            return Ok(BarcodePreprocessing {
                dir: dir.to_path_buf(),
                chunks: index.chunks,
                stats: saved.stats,
                total,
                orientation,
            });
        }
    }
    fs::create_dir_all(dir)?;
    let started = Instant::now();
    let orientation = if c.barcode.orientation == "auto" {
        detect_orientation(c)?
    } else if c.barcode.orientation == "rc" {
        "rc"
    } else {
        "forward"
    };
    let matcher = Matcher::load(c, orientation)?;
    eprintln!("event=barcode_start chunks={}", index.chunks);
    let starts = index.chunk_starts();
    let mut chunk = 0usize;
    let mut stats = vec![BarcodeStats::default(); index.chunks];
    let mut writer =
        BufWriter::with_capacity(1024 * 1024, File::create(dir.join("chunk-0000.tsv"))?);
    let barcode_read = c
        .barcode_read
        .as_deref()
        .context("barcode_read is required for the droplet barcode preprocessing")?;
    let decode_threads = c.barcode.matcher_threads.clamp(1, 4);
    let mut reader = open_barcode_read(barcode_read, decode_threads)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(c.barcode.matcher_threads as usize)
        .build()
        .context("build barcode matcher thread pool")?;
    let mut data: Vec<u8> = Vec::with_capacity(PREPASS_BLOCK * 160);
    let mut spans: Vec<(u32, u32, u32, u32)> = Vec::with_capacity(PREPASS_BLOCK);
    let (mut h, mut s, mut p, mut q) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut processed: u64 = 0;
    let mut next_report: u64 = 10_000_000;
    let mut eof = false;
    while !eof {
        while spans.len() < PREPASS_BLOCK {
            h.clear();
            s.clear();
            p.clear();
            q.clear();
            if reader.read_until(b'\n', &mut h)? == 0 {
                eof = true;
                break;
            }
            if reader.read_until(b'\n', &mut s)? == 0
                || reader.read_until(b'\n', &mut p)? == 0
                || reader.read_until(b'\n', &mut q)? == 0
            {
                bail!("truncated R2 FASTQ record")
            }
            if !h.starts_with(b"@") || !p.starts_with(b"+") {
                bail!("invalid R2 FASTQ record")
            }
            trim_newline(&mut h);
            trim_newline(&mut s);
            let name = name_of(&h);
            if chunk < starts.len() && name == starts[chunk].as_bytes() {
                flush_block(&pool, &matcher, &mut data, &mut spans, &mut stats[chunk], &mut writer)?;
                chunk += 1;
                writer.flush()?;
                writer = BufWriter::with_capacity(
                    1024 * 1024,
                    File::create(dir.join(format!("chunk-{chunk:04}.tsv")))?,
                );
            }
            // FASTQ headers start with '@'; SAM QNAMEs do not.
            let qname = name.strip_prefix(b"@").unwrap_or(name);
            let ns = data.len() as u32;
            data.extend_from_slice(qname);
            let ne = data.len() as u32;
            let ss = data.len() as u32;
            data.extend_from_slice(&s);
            let se = data.len() as u32;
            spans.push((ns, ne, ss, se));
        }
        let n = flush_block(&pool, &matcher, &mut data, &mut spans, &mut stats[chunk], &mut writer)?;
        processed += n as u64;
        if processed >= next_report {
            let seconds = started.elapsed().as_secs_f64().max(0.001);
            eprintln!(
                "event=barcode_progress reads={} elapsed_s={:.1} rate={:.0}/s",
                processed,
                seconds,
                processed as f64 / seconds
            );
            next_report = (processed / 10_000_000 + 1) * 10_000_000;
        }
    }
    writer.flush()?;
    if chunk + 1 != index.chunks {
        bail!(
            "R2 boundary read IDs matched {}/{} chunks; FASTQ files may be out of sync",
            chunk + 1,
            index.chunks
        )
    }
    let mut total = BarcodeStats::default();
    for s in &stats {
        total.add(s);
    }
    #[derive(Serialize)]
    struct Saved<'a> {
        orientation: &'a str,
        stats: &'a [BarcodeStats],
        slice_start: usize,
        slice_len: usize,
        max_mismatch: u8,
        matcher_version: u32,
    }
    fs::write(
        &stats_path,
        serde_json::to_vec_pretty(&Saved {
            orientation,
            stats: &stats,
            slice_start: c.barcode.slice_start,
            slice_len: c.barcode.slice_len,
            max_mismatch: c.barcode.max_mismatch,
            matcher_version: MATCHER_VERSION,
        })?,
    )?;
    eprintln!(
        "event=barcode_done total={} exact_match={} match_within_threshold={} ambiguous={} unassigned={} elapsed_seconds={:.3}",
        total.total,
        total.exact,
        total.corrected,
        total.ambiguous,
        total.unassigned + total.poly_g,
        started.elapsed().as_secs_f64()
    );
    Ok(BarcodePreprocessing {
        dir: dir.to_path_buf(),
        chunks: index.chunks,
        stats,
        total,
        orientation,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(s: &[u8]) -> [u8; 16] {
        let mut a = [0u8; 16];
        a[..s.len()].copy_from_slice(s);
        a
    }

    fn matcher() -> Matcher {
        Matcher::from_entries(
            vec![
                key(b"AAAACCCCGGGGTTTT"),
                key(b"CAAACCCCGGGGTTTT"),
            ],
            0,
            16,
            2,
            false,
        )
    }

    #[test]
    fn exact_unique_and_ambiguous_barcodes() {
        let m = matcher();
        assert_eq!(
            m.classify(b"AAAACCCCGGGGTTTT"),
            Match::Exact(key(b"AAAACCCCGGGGTTTT"))
        );
        assert_eq!(
            m.classify(b"AAAACCCCGGGGTTTA"),
            Match::Corrected(key(b"AAAACCCCGGGGTTTT"))
        );
        assert_eq!(m.classify(b"GAAACCCCGGGGTTTT"), Match::Ambiguous);
        assert_eq!(m.classify(b"TTTTTTTTTTTTTTTT"), Match::Unassigned);
    }

    #[test]
    fn n_base_is_a_mismatch_not_poly_g() {
        let m = matcher();
        assert_eq!(
            m.classify(b"AAAACCCCGGGGTTTN"),
            Match::Corrected(key(b"AAAACCCCGGGGTTTT"))
        );
        // Two Ns at the 3' end still leave the read within the two-mismatch
        // budget, so it is correctable.
        assert_eq!(
            m.classify(b"AAAACCCCGGGGTTNN"),
            Match::Corrected(key(b"AAAACCCCGGGGTTTT"))
        );
        // A read with no informative bases is a poly-G artifact.
        assert_eq!(m.classify(b"GGGGGGGGGGGGGGGG"), Match::PolyG);
    }

    #[test]
    fn two_mismatch_correction() {
        let m = matcher();
        assert_eq!(
            m.classify(b"AAAACCCCGGGGTTAA"),
            Match::Corrected(key(b"AAAACCCCGGGGTTTT"))
        );
    }

    #[test]
    fn rc_orientation_is_applied() {
        let bc = key(b"AAAACCCCGGGGTTTA");
        let m = Matcher::from_entries(vec![bc], 0, 16, 0, true);
        // Reverse complement of the whitelist barcode is the sequenced form.
        assert_eq!(m.classify(b"TAAACCCCGGGGTTTT"), Match::Exact(bc));
        // The whitelist orientation itself must not match under rc.
        assert_eq!(m.classify(b"AAAACCCCGGGGTTTA"), Match::Unassigned);
    }

    #[test]
    fn reader_handles_plain_and_multimember_gzip() {
        use std::io::{Read as _, Write as _};
        let dir = std::env::temp_dir().join(format!("rh_barcode_reader_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let records = [
            b"@r1\nAAAACCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n".as_slice(),
            b"@r2\nTTTTGGGGCCCCAAAA\n+\nIIIIIIIIIIIIIIII\n".as_slice(),
        ];

        let plain = dir.join("r2.fastq");
        let mut f = std::fs::File::create(&plain).unwrap();
        for r in records {
            f.write_all(r).unwrap();
        }
        drop(f);
        let mut text = String::new();
        open_barcode_read(&plain, 2)
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text.lines().count(), 8);

        let gz = dir.join("r2.fastq.gz");
        let mut out = std::fs::File::create(&gz).unwrap();
        for r in records {
            let mut enc =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            enc.write_all(r).unwrap();
            out.write_all(&enc.finish().unwrap()).unwrap();
        }
        drop(out);
        let mut text = String::new();
        open_barcode_read(&gz, 2)
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text.lines().count(), 8);
        std::fs::remove_dir_all(&dir).ok();
    }
}

