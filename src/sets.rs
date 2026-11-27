//! Which combinations of variants become output sequences.

use crate::plan::Plan;
use clap::ValueEnum;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Mode {
    /// Ladder: unmodified input, v1, v1+v2, v1+v2+v3, ... (spec order; parents always precede children).
    Cumulative,
    /// Each variant on its own (plus its ancestors, since a nested variant needs its parent).
    Single,
}

#[derive(Clone, Debug)]
pub struct OutSet {
    pub name: String,
    pub active: Vec<bool>,
}

/// `product_0000` is the unmodified input; the rest are `product_0001`, `product_0002`, ... in output order.
/// Which variants each product holds is recorded in `products.tsv`.
fn product_name(rank: usize, width: usize) -> String {
    format!("product_{rank:0width$}")
}

pub fn build_sets(plan: &Plan, mode: Mode, drop_original: bool) -> Vec<OutSet> {
    let n = plan.nodes.len();
    let mut sets: Vec<Vec<bool>> = Vec::new();
    if !drop_original {
        sets.push(vec![false; n]);
    }
    match mode {
        Mode::Cumulative => {
            for k in 1..=n {
                let mut a = vec![false; n];
                a[..k].iter_mut().for_each(|x| *x = true);
                sets.push(a);
            }
        }
        Mode::Single => {
            for i in 0..n {
                let mut a = vec![false; n];
                for j in plan.with_ancestors(i) {
                    a[j] = true;
                }
                sets.push(a);
            }
        }
    }
    let mut seen = std::collections::HashSet::new();
    let sets: Vec<Vec<bool>> = sets.into_iter().filter(|a| seen.insert(a.clone())).collect();
    // Padding fits the total including `product_0000`, whether or not it is written (`--drop_input`), so the same
    // spec gives the same file names either way and a rerun overwrites instead of leaving duplicates. Min 4 digits.
    let width = (sets.iter().filter(|a| a.iter().any(|&x| x)).count() + 1).to_string().len().max(4);
    let mut rank = 0;
    sets.into_iter()
        .map(|active| {
            if active.iter().any(|&x| x) {
                rank += 1;
            }
            OutSet { name: product_name(rank, width), active }
        })
        .collect()
}
