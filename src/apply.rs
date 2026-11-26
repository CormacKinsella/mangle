//! Realising a subset of the plan as sequence.
//!
//! A node is realised recursively: take its input region, replace each active child's span with that
//! child's realised output, then apply the node's own operation to the result. So a SNP inside an
//! inversion is applied first and then reverse-complemented with the rest of the inversion.
//!
//! Nested-semantics choices (documented in the README):
//! * child of `inv`: realised before the inversion, so it appears reverse-complemented;
//! * child of `dup`/`cnv`: lands in ONE extra copy (the last); the other copies stay pristine;
//! * child of `ins`: placed inside the generated inserted sequence;
//! * children of `snp`/`del` are rejected at parse time.

use crate::io::Reference;
use crate::plan::{Plan, Src};
use crate::spec::Kind;

pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' | b'u' => b'a',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'r' => b'y',
        b'y' => b'r',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        other => other, // N, S, W and anything unknown are self-complementary / left as-is
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

fn frame<'a>(plan: &'a Plan, refs: &'a Reference, n: usize) -> &'a [u8] {
    let node = &plan.nodes[n];
    match node.src {
        Src::Ref { contig, start } => &refs.seqs[contig][start..start + node.len],
        Src::Gen { node: g, off } => &plan.nodes[g].gen_seq[off..off + node.len],
    }
}

pub fn realise(plan: &Plan, refs: &Reference, active: &[bool], n: usize) -> Vec<u8> {
    let node = &plan.nodes[n];
    let f = frame(plan, refs, n);

    let mut kids: Vec<usize> = node.children.iter().copied().filter(|&c| active[c]).collect();
    kids.sort_by_key(|&c| (plan.nodes[c].off, plan.nodes[c].span()));
    let inner: Vec<u8> = if kids.is_empty() {
        f.to_vec()
    } else {
        let mut v = Vec::with_capacity(f.len());
        let mut cur = 0;
        for &c in &kids {
            let ch = &plan.nodes[c];
            v.extend_from_slice(&f[cur..ch.off]);
            v.extend_from_slice(&realise(plan, refs, active, c));
            cur = ch.off + ch.span();
        }
        v.extend_from_slice(&f[cur..]);
        v
    };

    match node.kind {
        Kind::Snp => vec![node.alt],
        Kind::Del => Vec::new(),
        Kind::Ins => inner,
        Kind::Inv => revcomp(&inner),
        Kind::Dup | Kind::Cnv => {
            let mut out = Vec::with_capacity(f.len() * node.copies);
            for _ in 0..node.copies - 1 {
                out.extend_from_slice(f);
            }
            out.extend_from_slice(&inner);
            out
        }
    }
}

/// Stream one contig of the output defined by `active`, without materialising the whole contig.
pub fn render_contig(plan: &Plan, refs: &Reference, active: &[bool], contig: usize, emit: &mut dyn FnMut(&[u8])) {
    let seq = &refs.seqs[contig];
    let mut top: Vec<usize> = (0..plan.nodes.len())
        .filter(|&i| active[i] && plan.nodes[i].parent.is_none() && plan.nodes[i].contig == contig)
        .collect();
    top.sort_by_key(|&i| plan.nodes[i].off);
    let mut cur = 0;
    for i in top {
        let n = &plan.nodes[i];
        emit(&seq[cur..n.off]);
        emit(&realise(plan, refs, active, i));
        cur = n.off + n.span();
    }
    emit(&seq[cur..]);
}

pub fn render_contig_vec(plan: &Plan, refs: &Reference, active: &[bool], contig: usize) -> Vec<u8> {
    let mut v = Vec::new();
    render_contig(plan, refs, active, contig, &mut |b| v.extend_from_slice(b));
    v
}
