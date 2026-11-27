//! Safety limits. Everything here runs before any file is written, so a refused run leaves no partial output.
//!
//! The output budget is computed from the plan and the list of products, not from the command line, so
//! new ways of choosing products (`--mode all`, a haplotype count, ploidy) are covered automatically: they
//! only change the list of products, or `haplotypes_per_product` below.

use crate::apply::contig_len;
use crate::io::Reference;
use crate::plan::Plan;
use crate::sets::OutSet;
use crate::spec::Kind;
use anyhow::{Result, bail};

/// Most copies a `cnv` may have. Fixed: a sanity limit, not a tuning knob.
pub const MAX_CNV_COPIES: usize = 500;
/// Most bases a single variant may add to the sequence (decimal GB). Fixed, because each variant is built
/// in memory; revisit if copies are ever streamed instead.
pub const MAX_VARIANT_GROWTH: u64 = 1_000_000_000;

/// Bases a variant adds to the sequence: `ins` and `dup` add their length, `cnv` adds `len x (copies-1)`.
pub fn variant_growth(kind: Kind, len: usize, copies: usize) -> u64 {
    match kind {
        Kind::Ins | Kind::Dup => len as u64,
        Kind::Cnv => (len as u64).saturating_mul(copies.saturating_sub(1) as u64),
        _ => 0,
    }
}

/// Reject a variant that is too big. `max_len` is the largest length it can take (the top of a range), so the
/// result never depends on the seed. `copies` is as given in the catalogue (`None` = default).
pub fn check_variant(kind: Kind, max_len: usize, copies: Option<usize>) -> Result<()> {
    let copies = match kind {
        Kind::Dup => 2,
        Kind::Cnv => copies.unwrap_or(3),
        _ => 1,
    };
    if kind == Kind::Cnv && copies > MAX_CNV_COPIES {
        bail!("cnv copies must be <= {MAX_CNV_COPIES} (got {copies})");
    }
    let growth = variant_growth(kind, max_len, copies);
    if growth > MAX_VARIANT_GROWTH {
        bail!(
            "{} would add up to {} to the sequence, over the per-variant limit of {}",
            kind.name(),
            human(growth),
            human(MAX_VARIANT_GROWTH)
        );
    }
    Ok(())
}

/// How the products will be written. Add fields here when output gets new dimensions.
#[derive(Clone, Copy, Debug)]
pub struct OutputShape {
    /// FASTA records written per contig for each product (1 now; the ploidy / haplotype count once supported).
    pub haplotypes_per_product: usize,
    /// `--combined` writes every record a second time into all.fa.
    pub combined: bool,
    pub line_width: usize,
}

/// Exact size in bytes of the FASTA output (headers and line breaks included) for `sets`.
/// Truth files (BED, VCF, catalogue) are small by comparison and not counted.
pub fn estimate_output_bytes(plan: &Plan, refs: &Reference, sets: &[OutSet], shape: &OutputShape) -> u64 {
    let width = shape.line_width.max(1) as u64;
    let mut total: u64 = 0;
    for s in sets {
        let mut per_haplotype: u64 = 0;
        for (c, name) in refs.names.iter().enumerate() {
            let contig = name.rsplit('#').next().unwrap_or(name);
            let len = contig_len(plan, refs, &s.active, c) as u64;
            let header = (s.name.len() + contig.len() + 5) as u64; // '>' + name + '#0#' + contig + '\n'
            per_haplotype = per_haplotype.saturating_add(header).saturating_add(len).saturating_add(len.div_ceil(width));
        }
        total = total.saturating_add(per_haplotype.saturating_mul(shape.haplotypes_per_product as u64));
    }
    if shape.combined { total.saturating_mul(2) } else { total }
}

/// Refuse to continue if the output would exceed `max_bytes`. Returns the estimate.
pub fn check_output_size(plan: &Plan, refs: &Reference, sets: &[OutSet], shape: &OutputShape, max_bytes: u64) -> Result<u64> {
    let est = estimate_output_bytes(plan, refs, sets, shape);
    if est > max_bytes {
        bail!(
            "output would be about {} ({} sequences x {} haplotype(s){}), over the limit of {}; \
             reduce the variants/products or raise --max_output_size",
            human(est),
            sets.len(),
            shape.haplotypes_per_product,
            if shape.combined { ", doubled by --combined" } else { "" },
            human(max_bytes)
        );
    }
    Ok(est)
}

pub fn human(b: u64) -> String {
    const U: [(&str, f64); 3] = [("GB", 1e9), ("MB", 1e6), ("kB", 1e3)];
    for (n, d) in U {
        if b as f64 >= d {
            return format!("{:.1} {n}", b as f64 / d);
        }
    }
    format!("{b} B")
}
