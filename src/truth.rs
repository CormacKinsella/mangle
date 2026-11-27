//! Truth files: the resolved catalogue table, a BED of variant locations, and an optional per-product VCF.
//!
//! All coordinates are in the ORIGINAL reference. For a variant nested inside another, POS is its position
//! in the parent's unmodified frame (so a SNP inside an inversion is reported where it sits on the
//! reference, not where it ends up after the reverse-complement). Variants inside inserted sequence have no
//! reference coordinate and are omitted from the VCF (see resolved_catalogue.tsv, column ins_offset).

use crate::apply::contig_len;
use crate::io::Reference;
use crate::plan::Plan;
use crate::sets::OutSet;
use crate::spec::Kind;
use anyhow::Result;
use std::io::Write;

pub fn write_resolved_catalogue(w: &mut impl Write, plan: &Plan, refs: &Reference) -> Result<()> {
    writeln!(w, "id\ttype\tcontig\tref_pos_1based\tlen\tparent\tcopies\tref\talt\tins_offset")?;
    for n in &plan.nodes {
        let parent = n.parent.map(|p| plan.nodes[p].id.as_str()).unwrap_or(".");
        let (r, a) = if n.kind == Kind::Snp {
            ((n.refb as char).to_string(), (n.alt as char).to_string())
        } else {
            (".".into(), ".".into())
        };
        let go = n.gen_off.map(|g| g.to_string()).unwrap_or_else(|| ".".into());
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            n.id,
            n.kind.name(),
            refs.names[n.contig],
            n.ref_pos + 1,
            n.len,
            parent,
            n.copies,
            r,
            a,
            go
        )?;
    }
    Ok(())
}

pub fn write_vcf(w: &mut impl Write, plan: &Plan, refs: &Reference, active: &[bool], set: &str) -> Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=mangle")?;
    writeln!(w, "##mangle_set={set}")?;
    for (n, s) in refs.names.iter().zip(&refs.seqs) {
        writeln!(w, "##contig=<ID={n},length={}>", s.len())?;
    }
    for k in ["DEL", "INS", "INV", "DUP", "CNV"] {
        writeln!(w, "##ALT=<ID={k},Description=\"{k}\">")?;
    }
    writeln!(w, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">")?;
    writeln!(w, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length (negative for deletions)\">")?;
    writeln!(w, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position in original reference\">")?;
    writeln!(w, "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Total copy number of region\">")?;
    writeln!(w, "##INFO=<ID=PARENT,Number=1,Type=String,Description=\"ID of enclosing variant\">")?;
    writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

    // Variants inside inserted sequence have no reference coordinate (their REF base would not exist in
    // the reference, which breaks `bcftools norm -f`), so they are left out here and kept in resolved_catalogue.tsv.
    // POS uses the VCF anchor convention for ins/del, so compute it first and sort on it.
    let mut rows: Vec<(usize, usize, usize, u8)> = Vec::new(); // (contig, pos, node, refb)
    for i in (0..plan.nodes.len()).filter(|&i| active[i] && plan.nodes[i].gen_off.is_none()) {
        let n = &plan.nodes[i];
        let seq = &refs.seqs[n.contig];
        let (pos, refb) = if matches!(n.kind, Kind::Ins | Kind::Del) {
            (n.ref_pos.max(1), seq[n.ref_pos.saturating_sub(1)])
        } else {
            (n.ref_pos + 1, seq[n.ref_pos])
        };
        rows.push((n.contig, pos, i, refb));
    }
    rows.sort();
    for (contig, pos, i, refb) in rows {
        let n = &plan.nodes[i];
        let mut info: Vec<String> = Vec::new();
        let alt: String = match n.kind {
            Kind::Snp => {
                info.push("SVTYPE=SNP".into());
                format!("{}\t{}", refb as char, n.alt as char)
            }
            k => {
                let t = k.name().to_ascii_uppercase();
                info.push(format!("SVTYPE={t}"));
                let svlen = if k == Kind::Del { -(n.len as i64) } else { n.len as i64 };
                info.push(format!("SVLEN={svlen}"));
                let end = match k {
                    Kind::Ins => pos,
                    Kind::Del => pos + n.len,
                    _ => pos + n.len - 1,
                };
                info.push(format!("END={end}"));
                if matches!(k, Kind::Dup | Kind::Cnv) {
                    info.push(format!("CN={}", n.copies));
                }
                format!("{}\t<{}>", refb as char, t)
            }
        };
        if let Some(p) = n.parent {
            info.push(format!("PARENT={}", plan.nodes[p].id));
        }
        writeln!(w, "{}\t{}\t{}\t{}\t.\t.\t{}", refs.names[contig], pos, n.id, alt, info.join(";"))?;
    }
    Ok(())
}

/// BED6 (+ type, parent) of every planned variant in ORIGINAL reference coordinates, 0-based half-open.
/// Insertions are zero-length (start == end). Variants inside inserted sequence have no reference
/// coordinate and are omitted. A nested variant is reported in its parent's unmodified frame.
pub fn write_bed(w: &mut impl Write, plan: &Plan, refs: &Reference) -> Result<()> {
    let mut rows: Vec<(usize, usize, usize)> = (0..plan.nodes.len())
        .filter(|&i| plan.nodes[i].gen_off.is_none())
        .map(|i| (plan.nodes[i].contig, plan.nodes[i].ref_pos, i))
        .collect();
    rows.sort();
    for (c, start, i) in rows {
        let n = &plan.nodes[i];
        let end = start + n.span();
        let parent = n.parent.map(|p| plan.nodes[p].id.as_str()).unwrap_or(".");
        writeln!(w, "{}\t{}\t{}\t{}\t0\t.\t{}\t{}", refs.names[c], start, end, n.id, n.kind.name(), parent)?;
    }
    Ok(())
}

/// Most variants per product for which `products.tsv` lists ids; if any product has more, every row uses ranges.
const MAX_LISTED_IDS: usize = 20;

/// 1-based catalogue positions as ranges: `[1,2,3,7,9,10]` -> `1-3,7,9-10`.
fn ranges(pos: &[usize]) -> String {
    let mut out: Vec<String> = Vec::new();
    let mut i = 0;
    while i < pos.len() {
        let mut j = i;
        while j + 1 < pos.len() && pos[j + 1] == pos[j] + 1 {
            j += 1;
        }
        out.push(if j == i { pos[i].to_string() } else { format!("{}-{}", pos[i], pos[j]) });
        i = j + 1;
    }
    out.join(",")
}

/// Which variants each product contains. Ids when every product has at most `MAX_LISTED_IDS` variants,
/// otherwise 1-based row numbers of `resolved_catalogue.tsv` as ranges (a cumulative ladder stays one range per row).
pub fn write_products(w: &mut impl Write, plan: &Plan, refs: &Reference, sets: &[OutSet]) -> Result<()> {
    let by_id = sets.iter().all(|s| s.active.iter().filter(|&&a| a).count() <= MAX_LISTED_IDS);
    writeln!(w, "product\tn_variants\tlength\t{}", if by_id { "variants" } else { "variant_rows" })?;
    for s in sets {
        let rows: Vec<usize> = s.active.iter().enumerate().filter(|(_, a)| **a).map(|(i, _)| i + 1).collect();
        let list = match (rows.is_empty(), by_id) {
            (true, _) => ".".to_string(),
            (false, true) => rows.iter().map(|&r| plan.nodes[r - 1].id.as_str()).collect::<Vec<_>>().join(","),
            (false, false) => ranges(&rows),
        };
        let len: usize = (0..refs.names.len()).map(|c| contig_len(plan, refs, &s.active, c)).sum();
        writeln!(w, "{}\t{}\t{}\t{}", s.name, rows.len(), len, list)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::ranges;

    #[test]
    fn ranges_compress_runs() {
        assert_eq!(ranges(&[1, 2, 3, 7, 9, 10]), "1-3,7,9-10");
        assert_eq!(ranges(&[5]), "5");
    }
}
