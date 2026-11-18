//! Supported pairs column layouts shared by every processing stage.
use anyhow::{bail, Context, Result};
use std::io::BufRead;

pub(crate) const BASE_COLUMNS: [&str; 8] = [
    "readID",
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
    "pair_type",
];
pub(crate) const FRAGMENT_COLUMNS: [&str; 6] = [
    "rfrag1",
    "rfrag_start1",
    "rfrag_end1",
    "rfrag2",
    "rfrag_start2",
    "rfrag_end2",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Format {
    Standard,
    Fragments,
    /// Base columns plus a trailing `CB` cell-barcode column.
    Cell,
}
impl Format {
    pub fn from_columns(line: &str) -> Result<Self> {
        let columns: Vec<_> = line
            .strip_prefix("#columns:")
            .context("missing #columns header")?
            .split_whitespace()
            .collect();
        if columns == BASE_COLUMNS {
            Ok(Self::Standard)
        } else if columns.len() == 14
            && columns[..8] == BASE_COLUMNS
            && columns[8..] == FRAGMENT_COLUMNS
        {
            Ok(Self::Fragments)
        } else if columns.len() == 9 && columns[..8] == BASE_COLUMNS && columns[8] == "CB" {
            Ok(Self::Cell)
        } else {
            bail!("unsupported pairs columns: expected standard 8 columns, 14 columns with restriction-fragment annotations, or 9 columns with CB")
        }
    }
    pub fn len(self) -> usize {
        match self {
            Self::Standard => 8,
            Self::Fragments => 14,
            Self::Cell => 9,
        }
    }
    pub fn columns_line(self) -> String {
        let mut s = format!("#columns: {}", BASE_COLUMNS.join(" "));
        if self == Self::Fragments {
            s.push(' ');
            s.push_str(&FRAGMENT_COLUMNS.join(" "));
        } else if self == Self::Cell {
            s.push_str(" CB");
        }
        s.push('\n');
        s
    }
    pub fn fields(self, line: &str) -> Result<Vec<&str>> {
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() != self.len() {
            bail!(
                "pairs row has {} fields but #columns declares {}",
                fields.len(),
                self.len()
            )
        }
        Ok(fields)
    }
}

/// Read a pairs header and leave the stream at the beginning of its body.
pub(crate) fn read_header<R: BufRead>(reader: &mut R) -> Result<(Vec<String>, Format)> {
    let mut header = Vec::new();
    let mut format = None;
    while reader.fill_buf()?.first() == Some(&b'#') {
        let mut line = String::new();
        reader.read_line(&mut line)?;
        if line.starts_with("#columns:") {
            if format.is_some() {
                bail!("duplicate #columns header")
            }
            format = Some(Format::from_columns(&line)?);
        }
        header.push(line);
    }
    Ok((header, format.context("missing #columns header")?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cell_layout_roundtrips_with_cb() {
        let line = "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type CB\n";
        let f = Format::from_columns(line).unwrap();
        assert_eq!(f, Format::Cell);
        assert_eq!(f.len(), 9);
        assert_eq!(f.columns_line(), line);
        let row = "r1\tchr1\t1\tchr1\t2\t+\t-\tUU\tACGTACGTACGTACGT";
        assert_eq!(f.fields(row).unwrap()[8], "ACGTACGTACGTACGT");
    }

    #[test]
    fn fragments_and_standard_are_unchanged() {
        let standard = "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n";
        assert_eq!(Format::from_columns(standard).unwrap(), Format::Standard);
        let fragments = "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type rfrag1 rfrag_start1 rfrag_end1 rfrag2 rfrag_start2 rfrag_end2\n";
        assert_eq!(Format::from_columns(fragments).unwrap(), Format::Fragments);
    }
}
