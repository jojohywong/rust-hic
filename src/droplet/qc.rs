//! Per-cell QC from the deduplicated pairs file.
use anyhow::{bail, Context, Result};
use serde::Serialize;
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::{BufRead, BufReader, Write},
    path::Path,
};

#[derive(Default, Clone, Serialize)]
pub(crate) struct CellStats {
    pub pairs: u64,
    pub cis: u64,
    pub trans: u64,
    pub cis_1kb: u64,
    pub cis_10kb: u64,
    pub uu: u64,
    pub ur: u64,
    pub ru: u64,
    pub other: u64,
}

#[derive(Serialize)]
pub(crate) struct Summary {
    pub total_pairs: u64,
    pub cells: usize,
    pub cis: u64,
    pub trans: u64,
    pub cis_1kb: u64,
    pub cis_10kb: u64,
    pub unassigned_pairs: u64,
}

pub(crate) fn per_cell(pairs: &Path, out_dir: &Path) -> Result<Summary> {
    fs::create_dir_all(out_dir)?;
    let mut cells: BTreeMap<String, CellStats> = BTreeMap::new();
    let mut total = CellStats::default();
    let mut unassigned = 0u64;
    let mut body = false;
    let mut columns: Option<usize> = None;
    for line in BufReader::with_capacity(1024 * 1024, File::open(pairs)?).lines() {
        let line = line.with_context(|| format!("read {}", pairs.display()))?;
        if !body {
            if line.starts_with("#columns:") {
                columns = Some(line.split_whitespace().count() - 1);
            }
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            body = true;
        }
        let f: Vec<_> = line.split('\t').collect();
        let n = columns.context("pairs file has no #columns header")?;
        if f.len() != n {
            bail!("pairs row has {} fields, expected {}", f.len(), n)
        }
        if n != 9 {
            bail!("per-cell QC requires the 9-column CB layout")
        }
        let cell = f[8];
        let p1: i64 = f[2].parse().context("pos1")?;
        let p2: i64 = f[4].parse().context("pos2")?;
        let st = if cell == "." {
            unassigned += 1;
            cells.entry("unassigned".into()).or_default()
        } else {
            cells.entry(cell.to_owned()).or_default()
        };
        st.pairs += 1;
        total.pairs += 1;
        if f[1] == f[3] {
            st.cis += 1;
            total.cis += 1;
            let d = (p2 - p1).abs();
            if d >= 1_000 {
                st.cis_1kb += 1;
                total.cis_1kb += 1;
            }
            if d >= 10_000 {
                st.cis_10kb += 1;
                total.cis_10kb += 1;
            }
        } else {
            st.trans += 1;
            total.trans += 1;
        }
        match f[7] {
            "UU" => {
                st.uu += 1;
                total.uu += 1;
            }
            "UR" => {
                st.ur += 1;
                total.ur += 1;
            }
            "RU" => {
                st.ru += 1;
                total.ru += 1;
            }
            _ => {
                st.other += 1;
                total.other += 1;
            }
        }
    }
    let mut rows: Vec<_> = cells.into_iter().collect();
    rows.sort_by(|a, b| b.1.pairs.cmp(&a.1.pairs).then_with(|| a.0.cmp(&b.0)));
    let mut out = fs::File::create(out_dir.join("per_cell.tsv"))?;
    writeln!(
        out,
        "cell_id\tpairs\tcis\ttrans\tcis_1kb\tcis_10kb\tUU\tUR\tRU\tother"
    )?;
    for (cell, s) in &rows {
        writeln!(
            out,
            "{cell}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            s.pairs, s.cis, s.trans, s.cis_1kb, s.cis_10kb, s.uu, s.ur, s.ru, s.other
        )?;
    }
    let summary = Summary {
        total_pairs: total.pairs,
        cells: rows.iter().filter(|(name, _)| name != "unassigned").count(),
        cis: total.cis,
        trans: total.trans,
        cis_1kb: total.cis_1kb,
        cis_10kb: total.cis_10kb,
        unassigned_pairs: unassigned,
    };
    fs::write(
        out_dir.join("summary.json"),
        serde_json::to_vec_pretty(&summary)?,
    )?;
    Ok(summary)
}
