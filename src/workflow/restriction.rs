//! Built-in restriction enzyme definitions, embedded from a versioned data file.
use std::sync::OnceLock;

#[derive(Clone, Debug)]
pub(crate) struct Enzyme {
    pub(crate) name: &'static str,
    pub(crate) site: &'static str,
    pub(crate) fst5: i32,
    pub(crate) fst3: i32,
    pub(crate) ovhg: i32,
    pub(crate) scd5: Option<i32>,
    pub(crate) scd3: Option<i32>,
}

fn enzymes() -> &'static [Enzyme] {
    static ENZYMES: OnceLock<Vec<Enzyme>> = OnceLock::new();
    ENZYMES.get_or_init(|| {
        include_str!("../../data/restriction-enzymes.tsv")
            .lines()
            .filter(|line| {
                !line.starts_with('#') && !line.starts_with("name\t") && !line.is_empty()
            })
            .map(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                assert_eq!(fields.len(), 7, "invalid embedded enzyme table row: {line}");
                let optional = |s: &str| {
                    if s == "." {
                        None
                    } else {
                        Some(s.parse().expect("invalid cut offset"))
                    }
                };
                Enzyme {
                    name: fields[0],
                    site: fields[1],
                    fst5: fields[2].parse().expect("invalid fst5"),
                    fst3: fields[3].parse().expect("invalid fst3"),
                    ovhg: fields[4].parse().expect("invalid overhang"),
                    scd5: optional(fields[5]),
                    scd3: optional(fields[6]),
                }
            })
            .collect()
    })
}

pub(crate) fn lookup(name: &str) -> Option<&'static Enzyme> {
    enzymes()
        .iter()
        .find(|enzyme| enzyme.name.eq_ignore_ascii_case(name))
}

/// Parse `Arima`, a single enzyme, or a `-`-separated cocktail.
pub(crate) fn parse_spec(spec: &str) -> Result<Vec<&'static Enzyme>, String> {
    let spec = spec.trim();
    if spec.is_empty() {
        return Err("digest enzyme must not be empty".into());
    }
    if spec.eq_ignore_ascii_case("arima") {
        return Ok(vec![lookup("MboI").unwrap(), lookup("HinfI").unwrap()]);
    }
    let mut result = Vec::new();
    for part in spec.split('-') {
        let part = part.trim();
        if part.is_empty() {
            return Err(format!("invalid digest enzyme list: {spec}"));
        }
        let enzyme = lookup(part)
            .ok_or_else(|| format!("unknown or unsupported restriction enzyme '{part}'"))?;
        if !result
            .iter()
            .any(|existing: &&Enzyme| existing.name.eq_ignore_ascii_case(enzyme.name))
        {
            result.push(enzyme);
        }
    }
    Ok(result)
}
