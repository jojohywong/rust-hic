//! Optional Tn5-shifted fragment table from primary BAM alignments.
//!
//! This is intended for ATAC-Hi-C / spatial libraries where the same BAM carries
//! Tn5 fragments as well as Hi-C contacts. It is off by default for Droplet Hi-C.
use anyhow::{bail, Context, Result};
use flate2::{write::GzEncoder, Compression};
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::Cigar;
use noodles_sam::alignment::record::data::field::Value;
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

const CB_TAG: [u8; 2] = [b'C', b'B'];

pub(crate) struct FragmentConfig {
    pub t5_shift: u64,
    pub min_length: u64,
    pub max_length: u64,
}

pub(crate) fn write(bam_path: &Path, output: &Path, cfg: &FragmentConfig) -> Result<u64> {
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut reader = bam::io::Reader::new(BufReader::with_capacity(
        1024 * 1024,
        File::open(bam_path).with_context(|| format!("open BAM {}", bam_path.display()))?,
    ));
    let header = reader.read_header().context("read BAM header")?;
    let refs: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| n.to_string())
        .collect();

    // (chrom index, start, end, barcode) -> count
    let mut fragments: BTreeMap<(u16, u32, u32, String), u32> = BTreeMap::new();
    for result in reader.records() {
        let record = result.context("read BAM record")?;
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
        {
            continue;
        }
        let Some(rid) = record
            .reference_sequence_id()
            .transpose()
            .context("read reference id")?
        else {
            continue;
        };
        let Some(start) = record
            .alignment_start()
            .transpose()
            .context("read alignment start")?
            .map(usize::from)
        else {
            continue;
        };
        let span = record.cigar().alignment_span().context("read CIGAR span")?;
        if span == 0 {
            continue;
        }
        let Some(cb) = record
            .data()
            .get(&CB_TAG)
            .transpose()
            .context("read CB tag")?
            .and_then(|v| match v {
                Value::String(s) => Some(s.to_string()),
                _ => None,
            })
        else {
            continue;
        };
        if cb.is_empty() || cb == "." {
            continue;
        }
        // BAM start is 1-based; fragments.tsv is 0-based, end-exclusive.
        let start0 = start.saturating_sub(1) as u64 + cfg.t5_shift;
        let end0 = (start as u64 - 1 + span as u64).saturating_sub(5);
        if end0 <= start0 {
            continue;
        }
        let len = end0 - start0;
        if len < cfg.min_length || len > cfg.max_length {
            continue;
        }
        if start0 > u32::MAX as u64 || end0 > u32::MAX as u64 || rid > u16::MAX as usize {
            bail!("fragment coordinates exceed the supported range")
        }
        *fragments
            .entry((rid as u16, start0 as u32, end0 as u32, cb))
            .or_insert(0) += 1;
    }

    let encoder = GzEncoder::new(
        BufWriter::with_capacity(1024 * 1024, File::create(output)?),
        Compression::fast(),
    );
    let mut out = BufWriter::new(encoder);
    writeln!(out, "chrom\tstart\tend\tbarcode\tcount")?;
    for ((rid, start, end, cb), count) in &fragments {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            refs.get(*rid as usize)
                .context("reference id out of range")?,
            start,
            end,
            cb,
            count
        )?;
    }
    out.flush()?;
    let encoder = out
        .into_inner()
        .map_err(|e| anyhow::anyhow!("flush fragment writer: {e}"))?;
    encoder.finish().context("finish fragments.tsv.gz")?;
    Ok(fragments.len() as u64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Read};

    #[test]
    fn writes_tn5_shifted_fragments_and_aggregates_duplicates() {
        let seq = "A".repeat(100);
        let qual = "I".repeat(100);
        let sam = format!(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n\
             r1\t65\tchr1\t100\t60\t100M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:AAAACCCCGGGGTTTT\n\
             r1\t129\tchr1\t100\t60\t100M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:AAAACCCCGGGGTTTT\n\
             r2\t65\tchr1\t100\t60\t100M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:AAAACCCCGGGGTTTT\n\
             r3\t2113\tchr1\t500\t60\t100M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:AAAACCCCGGGGTTTT\n"
        );
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("in.bam");
        crate::mapping::bam::sam_to_bam(Cursor::new(sam), &bam).unwrap();
        let out = dir.path().join("fragments.tsv.gz");
        let n = write(
            &bam,
            &out,
            &FragmentConfig {
                t5_shift: 4,
                min_length: 1,
                max_length: 5000,
            },
        )
        .unwrap();
        assert_eq!(n, 1);
        let mut text = String::new();
        flate2::read::GzDecoder::new(File::open(&out).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        let rows: Vec<_> = text.lines().filter(|l| !l.starts_with("chrom")).collect();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0], "chr1\t103\t194\tAAAACCCCGGGGTTTT\t3");
    }
}
