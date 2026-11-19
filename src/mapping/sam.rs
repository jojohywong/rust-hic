//! BWA SAM to selected Hi-C pairs (pairtools parse's 5unique policy).
use anyhow::{bail, Context, Result};
use std::io::{BufRead, Write};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Alignment {
    chrom: String,
    pos: i64,
    pos3: i64,
    strand: char,
    mapped: bool,
    unique: bool,
    kind: char,
    dist5: i64,
    read_span: i64,
    query_start: i64,
}

/// Convert bwa-mem3's structured fields using the same 5unique semantics as SAM.
pub(crate) fn from_fields(
    fields: &bwa_mem3_rs::AlignedFields<'_>,
    index: &bwa_mem3_rs::BwaIndex,
    min_mapq: u8,
) -> Result<Alignment> {
    let chrom = if fields.flag & 4 == 0 {
        if fields.tid < 0 || fields.tid as usize >= index.n_contigs() {
            bail!("mapped BWA record lacks reference coordinates")
        }
        index.contig_name(fields.tid as usize)
    } else {
        "*"
    };
    from_fields_with_chrom(fields, chrom, min_mapq)
}

fn from_fields_with_chrom(
    fields: &bwa_mem3_rs::AlignedFields<'_>,
    chrom: &str,
    min_mapq: u8,
) -> Result<Alignment> {
    let flag = fields.flag;
    if flag & 1 == 0 {
        bail!("unpaired BWA alignment")
    }
    let mapped = flag & 4 == 0;
    let unique = fields.mapq >= min_mapq;
    let mut a = Alignment::empty();
    a.mapped = mapped;
    a.unique = unique;
    a.kind = if !mapped {
        'N'
    } else if !unique {
        'M'
    } else {
        'U'
    };
    if mapped {
        if fields.tid < 0 || fields.pos < 0 || chrom == "*" {
            bail!("mapped BWA record lacks reference coordinates")
        }
        let (mut matched, mut ref_span, mut read_span, mut leading, mut trailing, mut soft_start) =
            (0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
        let mut before_alignment = true;
        for &c in fields.cigar {
            let n = i64::from(c >> 4);
            if n == 0 {
                bail!("invalid CIGAR")
            }
            match c & 15 {
                0 => {
                    matched += n;
                    ref_span += n;
                    read_span += n;
                    before_alignment = false;
                }
                1 => {
                    read_span += n;
                    before_alignment = false;
                }
                2 => {
                    ref_span += n;
                    before_alignment = false;
                }
                op @ (4 | 5) => {
                    if matched == 0 {
                        leading = n
                    } else {
                        trailing = n
                    }
                    if op == 4 && before_alignment {
                        soft_start += n
                    }
                }
                _ => bail!("unsupported BWA CIGAR operation"),
            }
        }
        if ref_span == 0 {
            bail!("invalid or empty mapped CIGAR")
        }
        let reverse = flag & 16 != 0;
        a.dist5 = if reverse { trailing } else { leading };
        a.read_span = read_span;
        a.query_start = soft_start;
        if unique {
            a.chrom = chrom.to_owned();
            let pos = i64::from(fields.pos) + 1;
            let end = pos + ref_span - 1;
            if reverse {
                a.pos = end;
                a.pos3 = pos;
            } else {
                a.pos = pos;
                a.pos3 = end;
            }
            a.strand = if reverse { '-' } else { '+' };
        }
    }
    Ok(a)
}
impl Alignment {
    fn empty() -> Self {
        Self {
            chrom: "!".into(),
            pos: 0,
            pos3: 0,
            strand: '-',
            mapped: false,
            unique: false,
            kind: 'N',
            dist5: 0,
            read_span: 0,
            query_start: 0,
        }
    }
}
fn parse_alignment(f: &[&str], min_mapq: u8) -> Result<(bool, Alignment)> {
    if f.len() < 11 {
        bail!("SAM record has fewer than 11 fields")
    }
    let flag: u16 = f[1].parse().context("invalid SAM flag")?;
    if flag & 0x1 == 0 {
        bail!("unpaired SAM record")
    }
    let mapq: u8 = f[4].parse().context("invalid SAM MAPQ")?;
    let mapped = flag & 4 == 0;
    let unique = mapq >= min_mapq;
    let mut a = Alignment::empty();
    a.mapped = mapped;
    a.unique = unique;
    a.kind = if !mapped {
        'N'
    } else if !unique {
        'M'
    } else {
        'U'
    };
    if mapped {
        let pos: i64 = f[3].parse().context("invalid SAM position")?;
        if pos < 1 || f[2] == "*" || f[5] == "*" {
            bail!("mapped SAM record lacks coordinates or CIGAR")
        }
        let (mut n, mut matched, mut ref_span, mut read_span, mut leading, mut trailing) =
            (0i64, 0, 0, 0, 0, 0);
        let mut soft_start = 0;
        let mut before_alignment = true;
        for ch in f[5].chars() {
            if let Some(d) = ch.to_digit(10) {
                n = n
                    .checked_mul(10)
                    .and_then(|x| x.checked_add(d as i64))
                    .context("CIGAR overflow")?;
                continue;
            }
            if n == 0 {
                bail!("invalid CIGAR")
            }
            match ch {
                'M' => {
                    matched += n;
                    ref_span += n;
                    read_span += n;
                    before_alignment = false;
                }
                'I' => {
                    read_span += n;
                    before_alignment = false;
                }
                'D' => {
                    ref_span += n;
                    before_alignment = false;
                }
                'S' | 'H' => {
                    if matched == 0 {
                        leading = n
                    } else {
                        trailing = n
                    }
                    if ch == 'S' && before_alignment {
                        soft_start += n
                    }
                }
                _ => bail!("unsupported CIGAR operation {ch}; expected BWA M/I/D/S/H"),
            }
            n = 0;
        }
        if n != 0 || ref_span == 0 {
            bail!("invalid or empty mapped CIGAR")
        }
        let reverse = flag & 16 != 0;
        a.dist5 = if reverse { trailing } else { leading };
        a.read_span = read_span;
        a.query_start = soft_start;
        if unique {
            a.chrom = f[2].to_owned();
            let end = pos
                .checked_add(ref_span - 1)
                .context("coordinate overflow")?;
            if reverse {
                a.pos = end;
                a.pos3 = pos;
            } else {
                a.pos = pos;
                a.pos3 = end;
            }
            a.strand = if reverse { '-' } else { '+' };
        }
    }
    Ok((flag & 64 != 0, a))
}
// Order alignments by their 5' distance and represent internal uncovered
// segments as gaps. Unaligned sequence before the first alignment is clipping,
// not an additional alignment in the walk.
fn normalize(mut a: Vec<Alignment>, max_gap: i64) -> Vec<Alignment> {
    if a.is_empty() {
        return vec![Alignment::empty()];
    }
    a.sort_by_key(|x| x.query_start);
    a.sort_by_key(|x| x.dist5);
    if a.len() == 1 && !a[0].mapped {
        return a;
    }
    let mut alignments = a.into_iter();
    let first = alignments.next().expect("nonempty alignment list");
    let mut normalized = Vec::new();
    // pairtools converts a 5' unaligned tail longer than max_gap into a null
    // alignment, which then participates in rescue (leading-clipping parity).
    if first.dist5 > max_gap {
        let mut gap = Alignment::empty();
        gap.dist5 = 0;
        gap.read_span = first.dist5;
        normalized.push(gap);
    }
    let mut end = first.dist5 + first.read_span;
    normalized.push(first);
    for alignment in alignments {
        if alignment.dist5 - end > max_gap {
            let mut gap = Alignment::empty();
            gap.dist5 = end;
            gap.read_span = alignment.dist5 - end;
            normalized.push(gap);
        }
        end = end.max(alignment.dist5 + alignment.read_span);
        normalized.push(alignment);
    }
    normalized
}
fn classify(
    mut a: Vec<Alignment>,
    mut b: Vec<Alignment>,
    max_molecule: i64,
    max_gap: i64,
) -> Option<(Alignment, Alignment)> {
    if a.is_empty() || b.is_empty() {
        return None;
    }
    a = normalize(a, max_gap);
    b = normalize(b, max_gap);
    let mut rescued = false;
    if (a.len() == 2 && b.len() == 1) || (a.len() == 1 && b.len() == 2) {
        let (chim, linear) = if a.len() == 2 {
            (&a, &b[0])
        } else {
            (&b, &a[0])
        };
        if linear.mapped && linear.unique {
            let c5 = &chim[0];
            let c3 = &chim[1];
            // pairtools rescue_walk: auto-rescue when the 3' chimeric alignment
            // is unmapped, or when the 5' chimeric alignment is not unique;
            // otherwise require same chromosome, opposite strands, inward
            // orientation and a molecule smaller than the size selection.
            let mut can = true;
            if c3.mapped && c5.unique {
                can = c3.chrom == linear.chrom && c3.strand != linear.strand;
                let delta = if linear.strand == '+' {
                    c3.pos - linear.pos
                } else {
                    linear.pos - c3.pos
                };
                can &= delta > 0 && delta + c3.dist5 + linear.dist5 <= max_molecule;
            }
            if can {
                if a.len() == 2 {
                    b[0].kind = 'R';
                } else {
                    a[0].kind = 'R';
                }
                rescued = true;
            }
        }
    }
    let pick = |v: Vec<Alignment>| -> Alignment {
        if v.len() == 1 {
            v[0].clone()
        } else {
            v.iter()
                .find(|x| x.mapped && x.unique)
                .unwrap_or(&v[0])
                .clone()
        }
    };
    // A rescued pair reports the 5' alignment of the chimeric side and the
    // rescued linear alignment, matching pairtools (MR/NR/UR/RU). Pairtools
    // orders the less informative alignment first, so R+M/R+N becomes M+R/N+R.
    let (mut x, mut y) = if rescued {
        (a[0].clone(), b[0].clone())
    } else {
        (pick(a), pick(b))
    };
    if x.kind == 'R' && matches!(y.kind, 'M' | 'N') {
        std::mem::swap(&mut x, &mut y);
    }
    // pairtools orders the less informative side first: N < M < U.
    let rank = |a: &Alignment| (u8::from(a.mapped), u8::from(a.unique));
    if rank(&x) > rank(&y) {
        std::mem::swap(&mut x, &mut y);
    }
    Some((x, y))
}
fn usable(a: &Alignment) -> bool {
    a.mapped && a.unique
}

fn pairs_overlap(
    left: (&Alignment, &Alignment),
    right: (&Alignment, &Alignment),
    mismatch: i64,
) -> bool {
    let (l5, l3) = left;
    let (r5, r3) = right;
    let left_overlap = match (usable(l5), usable(r3)) {
        (false, false) => true,
        (false, true) | (true, false) => false,
        (true, true) => l5.chrom == r3.chrom && l5.strand != r3.strand,
    };
    let right_overlap = match (usable(l3), usable(r5)) {
        (false, false) => true,
        (false, true) | (true, false) => false,
        (true, true) => l3.chrom == r5.chrom && l3.strand != r5.strand,
    };
    let same_pair = (l5.pos3 - r3.pos).abs() <= mismatch && (l3.pos - r5.pos3).abs() <= mismatch;
    left_overlap && right_overlap && same_pair
}

fn partial_overlap(a1: &Alignment, a2: &Alignment, max_insert: i64, mismatch: i64) -> bool {
    if !usable(a1) && !usable(a2) {
        return true;
    }
    if !usable(a1) || !usable(a2) {
        return false;
    }
    if a1.chrom != a2.chrom || a1.strand == a2.strand {
        return false;
    }
    let (min_size, distance) = if a1.strand == '+' {
        ((a1.pos3 - a1.pos).max(a2.pos - a2.pos3), a2.pos - a1.pos)
    } else {
        ((a1.pos - a1.pos3).max(a2.pos3 - a2.pos), a1.pos - a2.pos)
    };
    distance <= max_insert + mismatch && distance >= min_size - mismatch
}

fn flip_position(a: &mut Alignment) {
    std::mem::swap(&mut a.pos, &mut a.pos3);
}

fn flip_orientation(a: &mut Alignment) {
    a.strand = if a.strand == '+' { '-' } else { '+' };
}

fn format_pair(
    a1: &Alignment,
    a2: &Alignment,
    pair_type: &str,
    a2_pos3: Option<i64>,
) -> (Alignment, Alignment) {
    let mut x = a1.clone();
    let mut y = a2.clone();
    if let Some(pos) = a2_pos3 {
        y.pos3 = pos;
    }
    match pair_type {
        "R1" | "R1&2" => {
            flip_position(&mut y);
            flip_orientation(&mut y);
        }
        "R2" => {
            flip_position(&mut x);
            flip_orientation(&mut x);
        }
        _ => {}
    }
    (x, y)
}

/// Port of pairtools 1.1.3 `parse_complex_walk` with the settings used by
/// `--walks-policy all`: report_position="outer", report_orientation="pair".
fn complex_walk(
    algns1: Vec<Alignment>,
    algns2: Vec<Alignment>,
    max_insert: i64,
) -> Vec<(Alignment, Alignment)> {
    let n1 = algns1.len();
    let n2 = algns2.len();
    let mut out: Vec<(usize, Alignment, Alignment)> = Vec::new();
    let current_left_pair = 1usize;
    let mut current_right_pair = 1usize;
    let remaining_left_pairs = n1.saturating_sub(1);
    let mut remaining_right_pairs = n2.saturating_sub(1);
    let mut checked_right_pairs = 0usize;
    let mut is_overlap = false;

    if n1 >= 2 && n2 >= 2 {
        while remaining_left_pairs > checked_right_pairs && remaining_right_pairs > 0 {
            let pair1 = (
                &algns1[n1 - current_left_pair - 1],
                &algns1[n1 - current_left_pair],
            );
            let pair2 = (
                &algns2[n2 - current_right_pair - 1],
                &algns2[n2 - current_right_pair],
            );
            is_overlap = pairs_overlap(pair1, pair2, 3);
            if is_overlap {
                let mut last_left = current_left_pair;
                let mut last_right = current_right_pair;
                let mut checked_temp = checked_right_pairs;
                while is_overlap && checked_temp > 0 {
                    last_left += 1;
                    last_right -= 1;
                    let p1 = (&algns1[n1 - last_left - 1], &algns1[n1 - last_left]);
                    let p2 = (&algns2[n2 - last_right - 1], &algns2[n2 - last_right]);
                    is_overlap &= pairs_overlap(p1, p2, 3);
                    checked_temp -= 1;
                }
                if is_overlap {
                    current_right_pair += 1;
                    break;
                }
            }
            current_right_pair += 1;
            checked_right_pairs += 1;
            remaining_right_pairs -= 1;
        }
        if !is_overlap {
            current_right_pair = 1;
        }
    }

    let (mut last_left, mut last_right);
    if n2 == 0 {
        last_left = 1;
        last_right = 0;
    } else if current_right_pair == 1 {
        last_left = 1;
        last_right = 1;
        if partial_overlap(&algns1[n1 - 1], &algns2[n2 - 1], max_insert, 3) {
            if n1 >= 2 {
                let (x, y) = format_pair(&algns1[n1 - 2], &algns1[n1 - 1], "R1", None);
                out.push((n1 - 1, x, y));
                last_left = 2;
            }
            if n2 >= 2 {
                let (x, y) = format_pair(&algns2[n2 - 1], &algns2[n2 - 2], "R2", None);
                out.push((n1, x, y));
                last_right = 2;
            }
        } else {
            let (x, y) = format_pair(&algns1[n1 - 1], &algns2[n2 - 1], "R1-2", None);
            out.push((n1, x, y));
        }
    } else {
        last_left = current_right_pair;
        last_right = current_right_pair;
    }

    for i in 0..n1.saturating_sub(last_left) {
        let (x, y) = format_pair(&algns1[i], &algns1[i + 1], "R1", None);
        out.push((i + 1, x, y));
    }
    for offset in 0..current_right_pair.saturating_sub(1) {
        let idx_left = n1 - current_right_pair + offset;
        let idx_right = n2 - 1 - offset;
        if idx_left + 1 >= n1 || idx_right == 0 {
            continue;
        }
        let override_pos3 = Some(algns2[idx_right - 1].pos);
        let (x, y) = format_pair(
            &algns1[idx_left],
            &algns1[idx_left + 1],
            "R1&2",
            override_pos3,
        );
        out.push((idx_left + 1, x, y));
    }
    let order_end = current_right_pair.min(n2.saturating_sub(last_right));
    for i in 0..order_end {
        let shift = if current_right_pair > 1 { -1i64 } else { 0 };
        let index = (n1 as i64 + order_end as i64 - i as i64 + shift) as usize;
        let (x, y) = format_pair(&algns2[i + 1], &algns2[i], "R2", None);
        out.push((index, x, y));
    }
    out.sort_by_key(|(i, _, _)| *i);
    out.into_iter()
        .map(|(_, mut a, mut b)| {
            let rank = |x: &Alignment| (u8::from(x.mapped), u8::from(x.unique));
            if rank(&a) > rank(&b) {
                std::mem::swap(&mut a, &mut b);
            }
            (a, b)
        })
        .collect()
}

fn walk_pairs(
    left: Vec<Alignment>,
    right: Vec<Alignment>,
    max_insert: i64,
) -> Vec<(Alignment, Alignment)> {
    if left.len() <= 1 && right.len() <= 1 {
        let (mut a, mut b) = (left[0].clone(), right[0].clone());
        let rank = |x: &Alignment| (u8::from(x.mapped), u8::from(x.unique));
        if rank(&a) > rank(&b) {
            std::mem::swap(&mut a, &mut b);
        }
        return vec![(a, b)];
    }
    complex_walk(left, right, max_insert)
}

pub(crate) fn write_pair_all<W: Write>(
    name: &str,
    a: &mut Vec<Alignment>,
    b: &mut Vec<Alignment>,
    out: &mut W,
    counts: &mut Counts,
    max_molecule: i64,
    max_gap: i64,
) -> Result<()> {
    emit_all(name, a, b, out, counts, max_molecule, max_gap)
}

fn emit_all<W: Write>(
    name: &str,
    a: &mut Vec<Alignment>,
    b: &mut Vec<Alignment>,
    out: &mut W,
    counts: &mut Counts,
    max_molecule: i64,
    max_gap: i64,
) -> Result<()> {
    if name.is_empty() {
        return Ok(());
    }
    counts.input_pairs += 1;
    let left = normalize(std::mem::take(a), max_gap);
    let right = normalize(std::mem::take(b), max_gap);
    for (x, y) in walk_pairs(left, right, max_molecule) {
        writeln!(
            out,
            "{name}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}",
            x.chrom, x.pos, y.chrom, y.pos, x.strand, y.strand, x.kind, y.kind
        )?;
        if matches!((x.kind, y.kind), ('U', 'U') | ('U', 'R') | ('R', 'U')) {
            counts.selected_pairs += 1;
        }
    }
    Ok(())
}

#[derive(Default, Debug, serde::Serialize)]
pub struct Counts {
    pub input_pairs: u64,
    pub selected_pairs: u64,
}
fn emit<W: Write>(
    name: &str,
    a: &mut Vec<Alignment>,
    b: &mut Vec<Alignment>,
    out: &mut W,
    counts: &mut Counts,
    max_molecule: i64,
    max_gap: i64,
) -> Result<()> {
    if name.is_empty() {
        return Ok(());
    }
    counts.input_pairs += 1;
    if let Some((x, y)) = classify(std::mem::take(a), std::mem::take(b), max_molecule, max_gap) {
        writeln!(
            out,
            "{name}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}",
            x.chrom, x.pos, y.chrom, y.pos, x.strand, y.strand, x.kind, y.kind
        )?;
        if matches!((x.kind, y.kind), ('U', 'U') | ('U', 'R') | ('R', 'U')) {
            counts.selected_pairs += 1;
        }
    }
    Ok(())
}

pub(crate) fn write_pair<W: Write>(
    name: &str,
    a: &mut Vec<Alignment>,
    b: &mut Vec<Alignment>,
    out: &mut W,
    counts: &mut Counts,
    max_molecule: i64,
    max_gap: i64,
) -> Result<()> {
    emit(name, a, b, out, counts, max_molecule, max_gap)
}
pub fn convert<R: BufRead, W: Write>(
    input: R,
    mut out: W,
    min_mapq: u8,
    max_molecule: i64,
    max_gap: i64,
) -> Result<Counts> {
    if min_mapq == 0 || max_molecule < 0 || max_gap < 0 {
        bail!("invalid parsing thresholds")
    }
    writeln!(out,"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type")?;
    let (mut a, mut b) = (Vec::new(), Vec::new());
    let mut name = String::new();
    let mut counts = Counts::default();
    let mut body = false;
    for (i, line) in input.lines().enumerate() {
        let line = line.with_context(|| format!("read SAM line {}", i + 1))?;
        if line.starts_with('@') && !body {
            continue;
        }
        body = true;
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 11 || f[0].is_empty() || f[0] == "*" {
            bail!("invalid SAM record on line {}", i + 1)
        }
        if f[0] != name {
            emit(
                &name,
                &mut a,
                &mut b,
                &mut out,
                &mut counts,
                max_molecule,
                max_gap,
            )?;
            name.clear();
            name.push_str(f[0]);
        }
        let (first, aln) =
            parse_alignment(&f, min_mapq).with_context(|| format!("SAM line {}", i + 1))?;
        if first {
            a.push(aln)
        } else {
            b.push(aln)
        }
    }
    emit(
        &name,
        &mut a,
        &mut b,
        &mut out,
        &mut counts,
        max_molecule,
        max_gap,
    )?;
    out.flush()?;
    Ok(counts)
}

/// Multi-way (`all`) parse without CB, emitting one row per contact.
#[allow(dead_code)]
pub fn convert_all<R: BufRead, W: Write>(
    input: R,
    mut out: W,
    min_mapq: u8,
    max_molecule: i64,
    max_gap: i64,
) -> Result<Counts> {
    if min_mapq == 0 || max_molecule < 0 || max_gap < 0 {
        bail!("invalid parsing thresholds")
    }
    writeln!(out,"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type")?;
    let (mut a, mut b) = (Vec::new(), Vec::new());
    let mut name = String::new();
    let mut counts = Counts::default();
    let mut body = false;
    for (i, line) in input.lines().enumerate() {
        let line = line.with_context(|| format!("read SAM line {}", i + 1))?;
        if line.starts_with('@') && !body {
            continue;
        }
        body = true;
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 11 || f[0].is_empty() || f[0] == "*" {
            bail!("invalid SAM record on line {}", i + 1)
        }
        if f[0] != name {
            emit_all(
                &name,
                &mut a,
                &mut b,
                &mut out,
                &mut counts,
                max_molecule,
                max_gap,
            )?;
            name.clear();
            name.push_str(f[0]);
        }
        let (first, aln) =
            parse_alignment(&f, min_mapq).with_context(|| format!("SAM line {}", i + 1))?;
        if first {
            a.push(aln)
        } else {
            b.push(aln)
        }
    }
    emit_all(
        &name,
        &mut a,
        &mut b,
        &mut out,
        &mut counts,
        max_molecule,
        max_gap,
    )?;
    out.flush()?;
    Ok(counts)
}

#[cfg(test)]
mod structured_tests {
    use super::*;
    use bwa_mem3_rs::AlignedFields;

    fn sam_record(name: &str, flag: u16, pos: i64, mapq: u8, cigar: &str) -> String {
        format!("{name}\t{flag}\tchr1\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t*\t*\n")
    }

    #[test]
    fn clipped_single_alignments_remain_usable_on_either_mate_and_strand() {
        let mut sam = String::new();
        for (name, first_flag, first_cigar, second_flag, second_cigar) in [
            ("leading", 65, "21S79M", 129, "100M"),
            ("reverse", 81, "79M21S", 129, "100M"),
            ("mate", 65, "100M", 129, "30S70M"),
        ] {
            sam.push_str(&sam_record(name, first_flag, 100, 60, first_cigar));
            sam.push_str(&sam_record(name, second_flag, 1000, 60, second_cigar));
        }
        let mut output = Vec::new();
        let counts = convert(sam.as_bytes(), &mut output, 1, 750, 20).unwrap();
        let pairs = String::from_utf8(output).unwrap();
        assert_eq!(counts.input_pairs, 3);
        assert_eq!(counts.selected_pairs, 0);
        assert!(pairs.contains("leading\t!\t0\tchr1\t1000\t-\t+\tNR\n"));
        assert!(pairs.contains("reverse\t!\t0\tchr1\t1000\t-\t+\tNR\n"));
        assert!(pairs.contains("mate\t!\t0\tchr1\t100\t-\t+\tNR\n"));
    }

    #[test]
    fn leading_and_internal_gaps_are_materialized() {
        let alignment = |start, length| Alignment {
            chrom: "chr1".into(),
            pos: 100,
            pos3: 100 + length - 1,
            strand: '+',
            mapped: true,
            unique: true,
            kind: 'U',
            dist5: start,
            read_span: length,
            query_start: start,
        };
        let normalized = normalize(
            vec![alignment(90, 10), alignment(30, 10), alignment(60, 10)],
            20,
        );
        // Gaps equal to the threshold do not become separate alignments.
        assert_eq!(
            normalized
                .iter()
                .map(|a| (a.kind, a.dist5, a.read_span))
                .collect::<Vec<_>>(),
            [('N', 0, 30), ('U', 30, 10), ('U', 60, 10), ('U', 90, 10)]
        );
        let normalized = normalize(
            vec![alignment(0, 10), alignment(40, 10), alignment(80, 10)],
            20,
        );
        assert_eq!(
            normalized
                .iter()
                .map(|a| (a.kind, a.dist5, a.read_span))
                .collect::<Vec<_>>(),
            [
                ('U', 0, 10),
                ('N', 10, 30),
                ('U', 40, 10),
                ('N', 50, 30),
                ('U', 80, 10),
            ]
        );
    }

    #[test]
    fn unusable_first_chimeric_alignment_falls_back_to_unique_selection() {
        let mut sam = String::new();
        // The first alignment is non-unique, while the second is unique.
        // A rescue must not choose the non-unique first alignment.
        sam.push_str(&sam_record("fallback", 65, 100, 0, "40M60S"));
        sam.push_str(&sam_record("fallback", 2113, 300, 60, "40S60M"));
        sam.push_str(&sam_record("fallback", 129, 1000, 60, "100M"));
        // The same fallback applies when the chimeric read is mate 2.
        sam.push_str(&sam_record("mate-fallback", 65, 1000, 60, "100M"));
        sam.push_str(&sam_record("mate-fallback", 129, 100, 0, "40M60S"));
        sam.push_str(&sam_record("mate-fallback", 2177, 300, 60, "40S60M"));
        // A valid rescue must continue to produce U/R.
        sam.push_str(&sam_record("rescue", 65, 5000, 60, "40M60S"));
        sam.push_str(&sam_record("rescue", 2129, 1100, 60, "60M40S"));
        sam.push_str(&sam_record("rescue", 129, 1000, 60, "100M"));
        // A non-unique second alignment still falls back to 5unique.
        sam.push_str(&sam_record("second-nonunique", 65, 5000, 60, "40M60S"));
        sam.push_str(&sam_record("second-nonunique", 2113, 1100, 0, "40S60M"));
        sam.push_str(&sam_record("second-nonunique", 129, 1000, 60, "100M"));
        let mut output = Vec::new();
        let counts = convert(sam.as_bytes(), &mut output, 1, 750, 20).unwrap();
        let pairs = String::from_utf8(output).unwrap();
        assert_eq!(counts.input_pairs, 4);
        assert_eq!(counts.selected_pairs, 2);
        assert!(pairs.contains("fallback\t!\t0\tchr1\t1000\t-\t+\tMR\n"));
        assert!(pairs.contains("mate-fallback\t!\t0\tchr1\t1000\t-\t+\tMR\n"));
        assert!(pairs.contains("rescue\tchr1\t5000\tchr1\t1000\t+\t+\tUR\n"));
        assert!(pairs.contains("second-nonunique\tchr1\t5000\tchr1\t1000\t+\t+\tUU\n"));
    }

    #[test]
    fn structured_fields_match_sam_for_clips_indels_and_mapping_quality() {
        let cases: &[(u16, u8, &str, &[u32])] = &[
            (65, 60, "100M", &[(100 << 4)]),
            (65, 60, "21S79M", &[(21 << 4) | 4, (79 << 4)]),
            (81, 60, "79M21S", &[(79 << 4), (21 << 4) | 4]),
            (
                81,
                60,
                "5S45M2I3M4D45M5S",
                &[
                    (5 << 4) | 4,
                    (45 << 4),
                    (2 << 4) | 1,
                    (3 << 4),
                    (4 << 4) | 2,
                    (45 << 4),
                    (5 << 4) | 4,
                ],
            ),
            (65, 60, "10H90M", &[(10 << 4) | 5, (90 << 4)]),
            (81, 0, "100M", &[(100 << 4)]),
            (69, 0, "*", &[]),
        ];
        for &(flag, mapq, cigar, encoded) in cases {
            let fields = AlignedFields {
                tid: 0,
                pos: 99,
                next_tid: -1,
                next_pos: -1,
                tlen: 0,
                flag,
                bin: 0,
                mapq,
                query_start: 0,
                query_end: 100,
                cigar: encoded,
                nm: None,
                md: None,
                mc: None,
                mq: None,
                score: None,
                sub: None,
                rg: None,
                sa: None,
                pa: None,
                xa: None,
                hn: None,
                meth: None,
            };
            let flag_text = flag.to_string();
            let mapq_text = mapq.to_string();
            let sam = [
                "read", &flag_text, "chr1", "100", &mapq_text, cigar, "*", "0", "0", "*", "*",
            ];
            let (_, expected) = parse_alignment(&sam, 1).unwrap();
            assert_eq!(
                from_fields_with_chrom(&fields, "chr1", 1).unwrap(),
                expected,
                "{cigar} flag={flag}"
            );
        }
    }
}

#[cfg(test)]
mod all_policy_tests {
    use super::*;
    use std::io::Cursor;

    fn rec(name: &str, flag: u16, pos: i64, cigar: &str) -> String {
        let seq = "A".repeat(100);
        let qual = "I".repeat(100);
        format!(
            "{name}\t{flag}\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tCB:Z:AAAACCCCGGGGTTTT\n"
        )
    }

    #[test]
    fn nonchimeric_read_emits_one_pair() {
        let sam = format!(
            "@SQ\tSN:chr1\tLN:1000000\n{}{}",
            rec("r1", 65, 1000, "100M"),
            rec("r1", 129, 5000, "100M")
        );
        let mut out = Vec::new();
        let counts = convert_all(Cursor::new(sam), &mut out, 1, 750, 20).unwrap();
        assert_eq!(counts.input_pairs, 1);
        assert_eq!(counts.selected_pairs, 1);
        let body: Vec<_> = String::from_utf8(out)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(str::to_owned)
            .collect();
        assert_eq!(body.len(), 1);
        assert_eq!(body[0].split('\t').count(), 8);
    }

    #[test]
    fn unmapped_reads_emit_nn_like_pairtools() {
        let sam = format!(
            "@SQ\tSN:chr1\tLN:1000000\n{}{}",
            rec("r1", 77, 0, "*"),
            rec("r1", 141, 0, "*")
        );
        let mut out = Vec::new();
        let counts = convert_all(Cursor::new(sam), &mut out, 1, 750, 20).unwrap();
        assert_eq!(counts.input_pairs, 1);
        assert_eq!(counts.selected_pairs, 0);
        let body: Vec<_> = String::from_utf8(out)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(str::to_owned)
            .collect();
        assert!(body[0].ends_with("\tNN"), "{}", body[0]);
    }
}
