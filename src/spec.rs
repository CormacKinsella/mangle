//! Variant-catalogue parsing.
//!
//! Format: whitespace/tab separated, one variant per line, `#` comments, optional header.
//!
//! ```text
//! id    type       length   parent  [copies]
//! var1  inversion  1000     .
//! var2  snp        .        var1
//! var3  cnv        200-400  .       5
//! ```
//!
//! Parents must be declared before their children.

use anyhow::{Context, Result, bail};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Snp,
    Ins,
    Del,
    Inv,
    Dup,
    Cnv,
}

impl Kind {
    pub fn parse(s: &str) -> Option<Kind> {
        Some(match s.to_ascii_lowercase().as_str() {
            "snp" | "snv" => Kind::Snp,
            "ins" | "insertion" => Kind::Ins,
            "del" | "deletion" => Kind::Del,
            "inv" | "inversion" => Kind::Inv,
            "dup" | "duplication" => Kind::Dup,
            "cnv" => Kind::Cnv,
            _ => return None,
        })
    }

    pub fn name(&self) -> &'static str {
        match self {
            Kind::Snp => "snp",
            Kind::Ins => "ins",
            Kind::Del => "del",
            Kind::Inv => "inv",
            Kind::Dup => "dup",
            Kind::Cnv => "cnv",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Len {
    Fixed(usize),
    /// Inclusive range, drawn uniformly by the RNG.
    Range(usize, usize),
    Unspecified,
}

#[derive(Clone, Debug)]
pub struct VarSpec {
    pub id: String,
    pub kind: Kind,
    pub len: Len,
    pub parent: Option<String>,
    /// Total copy number of the region for `cnv` (default 3). `dup` is always 2.
    pub copies: Option<usize>,
}

/// Parse a size such as `500`, `20k`, `1_000` or `10g` (decimal suffixes k/m/g). Errors on overflow.
pub fn parse_size(s: &str) -> Result<usize> {
    let s = s.replace('_', "");
    let (num, mult) = match s.chars().last() {
        Some('k') | Some('K') => (&s[..s.len() - 1], 1_000usize),
        Some('m') | Some('M') => (&s[..s.len() - 1], 1_000_000),
        Some('g') | Some('G') => (&s[..s.len() - 1], 1_000_000_000),
        _ => (&s[..], 1),
    };
    let n: usize = num.parse().with_context(|| format!("bad size '{s}'"))?;
    n.checked_mul(mult).with_context(|| format!("size '{s}' is too large"))
}

fn parse_len(s: &str) -> Result<Len> {
    if matches!(s, "." | "-" | "length" | "NA") {
        return Ok(Len::Unspecified);
    }
    let range = s.split_once("..").or_else(|| s.split_once('-'));
    match range {
        Some((a, b)) if !a.is_empty() && !b.is_empty() => {
            let (a, b) = (parse_size(a)?, parse_size(b)?);
            if a > b {
                bail!("length range '{s}' has min > max");
            }
            Ok(Len::Range(a, b))
        }
        _ => Ok(Len::Fixed(parse_size(s)?)),
    }
}

/// Parse, discarding warnings (see [`parse_spec_warn`]).
pub fn parse_spec(text: &str) -> Result<Vec<VarSpec>> {
    Ok(parse_spec_warn(text)?.0)
}

/// Parse, also returning non-fatal warnings (e.g. a SNP given length 20 is used as 1 bp).
pub fn parse_spec_warn(text: &str) -> Result<(Vec<VarSpec>, Vec<String>)> {
    let mut warnings: Vec<String> = Vec::new();
    let mut out: Vec<VarSpec> = Vec::new();
    let mut seen_row = false;
    let mut ids: std::collections::HashSet<String> = std::collections::HashSet::new();
    for (ln, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if !seen_row {
            seen_row = true;
            if f[0].eq_ignore_ascii_case("id") || f[0].eq_ignore_ascii_case("identity") || f.get(1).is_some_and(|t| t.eq_ignore_ascii_case("type")) {
                continue; // header
            }
        }
        let at = || format!("spec line {}", ln + 1);
        if f.len() < 3 {
            bail!("{}: need at least id, type, length", at());
        }
        let id = f[0].to_string();
        if !id
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.'))
        {
            bail!("{}: id '{id}' may only contain [A-Za-z0-9_.-]", at());
        }
        if id.eq_ignore_ascii_case("original") || id.eq_ignore_ascii_case("combined") {
            bail!("{}: '{id}' is reserved (original = unmodified sequence, combined = the --combined output file)", at());
        }
        if out.len() >= crate::limits::MAX_VARIANTS {
            bail!("{}: more than {} variants in the spec", at(), crate::limits::MAX_VARIANTS);
        }
        if !ids.insert(id.to_ascii_lowercase()) {
            bail!("{}: duplicate id '{id}' (ids are compared ignoring case)", at());
        }
        let kind = Kind::parse(f[1]).with_context(|| format!("{}: unknown type '{}'", at(), f[1]))?;
        let mut len = parse_len(f[2]).with_context(at)?;
        let parent = match f.get(3).copied() {
            None | Some(".") | Some("-") | Some("0") | Some("none") => None,
            Some(p) => Some(p.to_string()),
        };
        let mut copies = match f.get(4).copied() {
            None | Some(".") => None,
            Some(c) => Some(c.parse::<usize>().with_context(|| format!("{}: bad copies '{c}'", at()))?),
        };

        match kind {
            Kind::Snp => {
                if !matches!(len, Len::Unspecified | Len::Fixed(1)) {
                    warnings.push(format!("{}: snp '{id}' is always 1 bp; ignoring length '{}'", at(), f[2]));
                }
                len = Len::Unspecified;
            }
            _ => {
                if len == Len::Unspecified {
                    bail!("{}: type {} needs a length", at(), kind.name());
                }
                if matches!(len, Len::Fixed(0) | Len::Range(0, _)) {
                    bail!("{}: length must be >= 1", at());
                }
            }
        }
        match (kind, copies) {
            (Kind::Cnv, Some(c)) if c < 2 => bail!("{}: cnv copies must be >= 2", at()),
            (Kind::Dup, Some(c)) if c != 2 => {
                warnings.push(format!("{}: dup '{id}' always has 2 copies; ignoring copies={c} (use cnv)", at()));
                copies = None;
            }
            (Kind::Cnv | Kind::Dup, _) => {}
            (k, Some(c)) => {
                warnings.push(format!("{}: copies={c} only applies to cnv; ignoring for {} '{id}'", at(), k.name()));
                copies = None;
            }
            _ => {}
        }
        let max_len = match len {
            Len::Fixed(n) | Len::Range(_, n) => n,
            Len::Unspecified => 1,
        };
        crate::limits::check_variant(kind, max_len, copies).with_context(|| format!("{}: variant '{id}'", at()))?;
        out.push(VarSpec { id, kind, len, parent, copies });
    }

    let kinds: HashMap<&str, Kind> = out.iter().map(|v| (v.id.as_str(), v.kind)).collect();
    let pos: HashMap<&str, usize> = out.iter().enumerate().map(|(i, v)| (v.id.as_str(), i)).collect();
    let mut depth: Vec<usize> = vec![1; out.len()]; // top-level = 1; parents precede children, so this is final
    for (i, v) in out.iter().enumerate() {
        if let Some(p) = &v.parent {
            let pi = *pos
                .get(p.as_str())
                .with_context(|| format!("variant '{}': unknown parent '{p}'", v.id))?;
            if pi >= i {
                bail!("variant '{}': parent '{p}' must be declared before it", v.id);
            }
            depth[i] = depth[pi] + 1;
            if depth[i] > crate::limits::MAX_NEST_DEPTH {
                bail!("variant '{}': nested {} levels deep, over the limit of {}", v.id, depth[i], crate::limits::MAX_NEST_DEPTH);
            }
            match kinds[p.as_str()] {
                Kind::Snp | Kind::Del => bail!(
                    "variant '{}': cannot nest inside {} '{p}'",
                    v.id,
                    kinds[p.as_str()].name()
                ),
                _ => {}
            }
        }
    }
    if out.is_empty() {
        bail!("spec contains no variants");
    }
    Ok((out, warnings))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_user_example() {
        let s = parse_spec("identity type length nest\nvar1 inversion 10 0\nvar2 snp length var1\n").unwrap();
        assert_eq!(s.len(), 2);
        assert_eq!(s[1].parent.as_deref(), Some("var1"));
        assert_eq!(s[1].len, Len::Unspecified);
    }

    #[test]
    fn ranges_and_suffixes() {
        assert_eq!(parse_len("20-50").unwrap(), Len::Range(20, 50));
        assert_eq!(parse_len("1k..2k").unwrap(), Len::Range(1000, 2000));
        assert_eq!(parse_len("5").unwrap(), Len::Fixed(5));
    }

    #[test]
    fn fixed_length_types_warn_and_continue() {
        let (v, w) = parse_spec_warn("a snp 20 .\nb dup 10 . 5\nc inv 10 . 3\n").unwrap();
        assert_eq!(v[0].len, Len::Unspecified);
        assert_eq!((v[1].copies, v[2].copies), (None, None));
        assert_eq!(w.len(), 3);
    }

    #[test]
    fn rejects_bad_nesting() {
        assert!(parse_spec("a snp . .\nb inv 5 a\n").is_err());
        assert!(parse_spec("a del 10 .\nb snp . a\n").is_err());
        assert!(parse_spec("b inv 5 a\na inv 10 .\n").is_err());
    }

    #[test]
    fn nesting_depth_is_capped() {
        let chain = |n: usize| {
            let mut t = String::from("v0 inv 1k .\n");
            for i in 1..n {
                t += &format!("v{i} inv 1k v{}\n", i - 1);
            }
            t
        };
        assert!(parse_spec(&chain(10)).is_ok());
        assert!(parse_spec(&chain(11)).is_err());
    }
}
