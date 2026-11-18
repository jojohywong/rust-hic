//! Explicit pair-type selection after deduplication.
//!
//! pairtools dedup clusters and retains all pair types. rust-hic's historical
//! output is restricted to UU/UR/RU; this stage makes that policy explicit,
//! configurable and accounted for in the report.
use super::format::read_header;
use anyhow::{bail, Result};
use serde::Serialize;
use std::{
    collections::BTreeMap,
    env,
    fs::{self, OpenOptions},
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
};

#[derive(Default, Serialize)]
struct Stats {
    input_records: u64,
    output_records: u64,
    removed_records: u64,
    removed_by_type: BTreeMap<String, u64>,
}

fn run() -> Result<()> {
    let args: Vec<_> = env::args().collect();
    if args.len() != 3 {
        bail!("usage: select_pairs TYPES STATS.json < dedup.pairs > selected.pairs")
    }
    let allowed: Vec<&str> = args[1].split(',').filter(|x| !x.is_empty()).collect();
    if allowed.is_empty() || allowed.iter().any(|t| t.len() != 2 || !t.is_ascii()) {
        bail!("invalid output pair types: {}", args[1])
    }
    let stats_path = Path::new(&args[2]);
    if stats_path.exists() {
        bail!("stats output exists")
    }
    let stdin = io::stdin();
    let mut input = BufReader::with_capacity(1024 * 1024, stdin.lock());
    let mut out = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    let (header, _) = read_header(&mut input)?;
    for h in &header {
        out.write_all(h.as_bytes())?;
    }
    let mut stats = Stats::default();
    let mut line = String::new();
    loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\r', '\n']);
        let f: Vec<_> = trimmed.split('\t').collect();
        if f.len() < 8 {
            bail!("pairs row has fewer than 8 columns")
        }
        stats.input_records += 1;
        if allowed.contains(&f[7]) {
            out.write_all(trimmed.as_bytes())?;
            out.write_all(b"\n")?;
            stats.output_records += 1;
        } else {
            stats.removed_records += 1;
            *stats.removed_by_type.entry(f[7].to_owned()).or_insert(0) += 1;
        }
    }
    out.flush()?;
    let partial = stats_path.with_extension(format!("partial.{}", std::process::id()));
    let mut f = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&partial)?;
    serde_json::to_writer_pretty(&mut f, &stats)?;
    writeln!(f)?;
    f.sync_all()?;
    fs::hard_link(&partial, stats_path)?;
    fs::remove_file(&partial)?;
    eprintln!(
        "event=select_done input_records={} output_records={} removed_records={} removed_types={}",
        stats.input_records,
        stats.output_records,
        stats.removed_records,
        stats.removed_by_type.len()
    );
    Ok(())
}

pub(crate) fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
