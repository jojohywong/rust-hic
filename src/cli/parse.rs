//! Hidden SAM→pairs parser used to compare against pairtools on identical input.
use super::ParseSamArgs;
use crate::mapping::sam;
use anyhow::{bail, Result};
use std::io::{self, BufWriter, Write};

pub(crate) fn run(args: ParseSamArgs) -> Result<()> {
    let stdin = io::stdin();
    let input = io::BufReader::with_capacity(1024 * 1024, stdin.lock());
    let stdout = io::stdout();
    let mut out = BufWriter::with_capacity(1024 * 1024, stdout.lock());
    match (args.walks_policy.as_str(), args.cb) {
        ("5unique", false) => sam::convert(
            input,
            &mut out,
            args.min_mapq,
            args.max_molecule_size,
            args.max_inter_align_gap,
        )?,
        ("5unique", true) => sam::convert_cb(
            input,
            &mut out,
            args.min_mapq,
            args.max_molecule_size,
            args.max_inter_align_gap,
        )?,
        ("all", false) => sam::convert_all(
            input,
            &mut out,
            args.min_mapq,
            args.max_molecule_size,
            args.max_inter_align_gap,
        )?,
        ("all", true) => sam::convert_all_cb(
            input,
            &mut out,
            args.min_mapq,
            args.max_molecule_size,
            args.max_inter_align_gap,
        )?,
        (other, _) => bail!("walks policy must be 5unique or all, got {other:?}"),
    };
    out.flush()?;
    Ok(())
}
