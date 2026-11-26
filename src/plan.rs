//! Resolving the random plan.
//!
//! Every position, length, allele and inserted sequence is drawn **once**, here, before any output is
//! written. Each output sequence is then just a subset of this single plan, so a variant sits at the same
//! place in every sequence that contains it.
//!
//! Coordinate frames: a top-level variant lives in reference coordinates. A nested variant lives at an
//! `off`set inside its parent's *input* region (the parent's region before the parent's own operation is
//! applied). Siblings never overlap, and children are placed independently of which of them end up active.

use crate::io::Reference;
use crate::spec::{Kind, Len, VarSpec};
use anyhow::{Result, bail};
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::collections::HashMap;

const MAX_TRIES: usize = 20_000;

/// Where the bytes of a node's input region come from.
#[derive(Clone, Copy, Debug)]
pub enum Src {
    Ref { contig: usize, start: usize },
    /// A slice of the generated sequence of insertion node `node`.
    Gen { node: usize, off: usize },
}

#[derive(Clone, Debug)]
pub struct Node {
    pub id: String,
    pub kind: Kind,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    /// Total copy number (dup = 2, cnv = n, otherwise 1).
    pub copies: usize,
    /// Length of the node's input region (for `ins`: length of the inserted sequence).
    pub len: usize,
    /// Offset in the parent's input region; for top-level nodes, the reference start.
    pub off: usize,
    pub contig: usize,
    pub src: Src,
    /// Inserted sequence (`ins` only).
    pub gen_seq: Vec<u8>,
    pub refb: u8,
    pub alt: u8,
    /// Reference coordinate (0-based) of the node, or of the nearest enclosing insertion point when the
    /// node sits inside inserted sequence.
    pub ref_pos: usize,
    /// Offset inside the nearest enclosing inserted sequence, if any.
    pub gen_off: Option<usize>,
}

impl Node {
    /// Bases of the parent's input region consumed by this node (insertions consume none).
    pub fn span(&self) -> usize {
        if self.kind == Kind::Ins { 0 } else { self.len }
    }
}

#[derive(Clone, Debug)]
pub struct Plan {
    pub nodes: Vec<Node>,
}

impl Plan {
    /// `i` plus all of its ancestors.
    pub fn with_ancestors(&self, mut i: usize) -> Vec<usize> {
        let mut v = vec![i];
        while let Some(p) = self.nodes[i].parent {
            v.push(p);
            i = p;
        }
        v
    }
}

fn pick_alt(refb: u8, rng: &mut ChaCha8Rng) -> u8 {
    let up = refb.to_ascii_uppercase();
    let c: Vec<u8> = b"ACGT".iter().copied().filter(|&b| b != up).collect();
    let a = c[rng.gen_range(0..c.len())];
    if refb.is_ascii_lowercase() { a.to_ascii_lowercase() } else { a }
}

fn base_at(nodes: &[Node], refs: &Reference, src: Src) -> u8 {
    match src {
        Src::Ref { contig, start } => refs.seqs[contig][start],
        Src::Gen { node, off } => nodes[node].gen_seq[off],
    }
}

pub fn resolve(specs: &[VarSpec], refs: &Reference, seed: u64, min_gap: usize) -> Result<Plan> {
    if refs.seqs.is_empty() {
        bail!("reference has no sequences");
    }
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let index: HashMap<&str, usize> = specs.iter().enumerate().map(|(i, s)| (s.id.as_str(), i)).collect();

    let mut nodes: Vec<Node> = Vec::with_capacity(specs.len());
    for (i, s) in specs.iter().enumerate() {
        let len = match (s.kind, &s.len) {
            (Kind::Snp, _) => 1,
            (_, Len::Fixed(n)) => *n,
            (_, Len::Range(a, b)) => rng.gen_range(*a..=*b),
            (_, Len::Unspecified) => bail!("variant '{}' needs a length", s.id),
        };
        let copies = match s.kind {
            Kind::Dup => 2,
            Kind::Cnv => s.copies.unwrap_or(3),
            _ => 1,
        };
        let gen_seq: Vec<u8> = if s.kind == Kind::Ins {
            (0..len).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect()
        } else {
            Vec::new()
        };
        nodes.push(Node {
            id: s.id.clone(),
            kind: s.kind,
            parent: s.parent.as_ref().map(|p| index[p.as_str()]),
            children: Vec::new(),
            copies,
            len,
            off: 0,
            contig: 0,
            src: if s.kind == Kind::Ins { Src::Gen { node: i, off: 0 } } else { Src::Ref { contig: 0, start: 0 } },
            gen_seq,
            refb: b'N',
            alt: b'N',
            ref_pos: 0,
            gen_off: None,
        });
    }
    for i in 0..nodes.len() {
        if let Some(p) = nodes[i].parent {
            nodes[p].children.push(i);
        }
    }

    // Pass 1: top-level variants, placed on the reference without overlap.
    let mut footprints: Vec<(usize, usize, usize)> = Vec::new();
    for i in 0..nodes.len() {
        if nodes[i].parent.is_none() {
            place_top(&mut nodes, i, refs, &mut rng, &mut footprints, min_gap)?;
        }
    }
    // Pass 2: children, one parent at a time. Parents precede children in the spec, so a parent's own
    // location is always known by the time its children are placed.
    for p in 0..nodes.len() {
        place_children(&mut nodes, p, refs, &mut rng)?;
    }
    Ok(Plan { nodes })
}

fn place_top(
    nodes: &mut [Node],
    i: usize,
    refs: &Reference,
    rng: &mut ChaCha8Rng,
    foot: &mut Vec<(usize, usize, usize)>,
    min_gap: usize,
) -> Result<()> {
    let fp = nodes[i].span().max(1);
    let min_start = if matches!(nodes[i].kind, Kind::Ins | Kind::Del) { 1 } else { 0 };
    let total: usize = refs.seqs.iter().map(|s| s.len()).sum();
    for _ in 0..MAX_TRIES {
        // contig weighted by length
        let mut r = rng.gen_range(0..total);
        let mut c = 0;
        while r >= refs.seqs[c].len() {
            r -= refs.seqs[c].len();
            c += 1;
        }
        let clen = refs.seqs[c].len();
        if clen < fp + min_start {
            continue;
        }
        let start = rng.gen_range(min_start..=clen - fp);
        if refs.seqs[c][start..start + fp].iter().any(|b| matches!(b, b'N' | b'n')) {
            continue;
        }
        if foot
            .iter()
            .any(|&(fc, fs, fe)| fc == c && start < fe + min_gap && fs < start + fp + min_gap)
        {
            continue;
        }
        foot.push((c, start, start + fp));
        let n = &mut nodes[i];
        n.contig = c;
        n.off = start;
        n.ref_pos = start;
        n.gen_off = None;
        if n.kind != Kind::Ins {
            n.src = Src::Ref { contig: c, start };
        }
        if n.kind == Kind::Snp {
            n.refb = refs.seqs[c][start];
            n.alt = pick_alt(n.refb, rng);
        }
        return Ok(());
    }
    bail!(
        "could not place variant '{}' (len {}) without overlap after {} tries: reference too small or too crowded (see --min_gap)",
        nodes[i].id,
        nodes[i].len,
        MAX_TRIES
    )
}

fn place_children(nodes: &mut Vec<Node>, p: usize, refs: &Reference, rng: &mut ChaCha8Rng) -> Result<()> {
    let kids = nodes[p].children.clone();
    if kids.is_empty() {
        return Ok(());
    }
    let consumed: usize = kids.iter().map(|&k| nodes[k].span()).sum();
    let plen = nodes[p].len;
    if consumed > plen {
        bail!(
            "children of '{}' need {} bp but '{}' is only {} bp (children: {})",
            nodes[p].id,
            consumed,
            nodes[p].id,
            plen,
            kids.iter().map(|&k| nodes[k].id.as_str()).collect::<Vec<_>>().join(", ")
        );
    }
    // Draw gaps instead of rejection sampling: always succeeds when the children fit.
    let leftover = plen - consumed;
    let mut order = kids;
    order.shuffle(rng);
    let mut cuts: Vec<usize> = (0..order.len()).map(|_| rng.gen_range(0..=leftover)).collect();
    cuts.sort_unstable();
    let mut used = 0;
    for (&k, cut) in order.iter().zip(cuts) {
        let off = cut + used;
        used += nodes[k].span();

        let (psrc, pkind, pref, pgen, pcontig) =
            (nodes[p].src, nodes[p].kind, nodes[p].ref_pos, nodes[p].gen_off, nodes[p].contig);
        let (ref_pos, gen_off) = if pkind == Kind::Ins {
            (pref, Some(off))
        } else if let Some(g) = pgen {
            (pref, Some(g + off))
        } else {
            (pref + off, None)
        };
        let src = if nodes[k].kind == Kind::Ins {
            Src::Gen { node: k, off: 0 }
        } else {
            match psrc {
                Src::Ref { contig, start } => Src::Ref { contig, start: start + off },
                Src::Gen { node, off: o } => Src::Gen { node, off: o + off },
            }
        };
        {
            let n = &mut nodes[k];
            n.off = off;
            n.contig = pcontig;
            n.ref_pos = ref_pos;
            n.gen_off = gen_off;
            n.src = src;
        }
        if nodes[k].kind == Kind::Snp {
            let b = base_at(nodes, refs, nodes[k].src);
            nodes[k].refb = b;
            nodes[k].alt = pick_alt(b, rng);
        }
    }
    Ok(())
}
