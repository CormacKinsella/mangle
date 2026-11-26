//! Which combinations of variants become output sequences.

use crate::plan::Plan;
use clap::ValueEnum;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Mode {
    /// Ladder: original, v1, v1+v2, v1+v2+v3, ... (spec order; parents always precede children).
    Cumulative,
    /// Each variant on its own (plus its ancestors, since a nested variant needs its parent).
    Single,
}

#[derive(Clone, Debug)]
pub struct OutSet {
    pub name: String,
    pub active: Vec<bool>,
}

fn name_of(plan: &Plan, active: &[bool]) -> String {
    let ids: Vec<&str> = plan
        .nodes
        .iter()
        .zip(active)
        .filter(|(_, a)| **a)
        .map(|(n, _)| n.id.as_str())
        .collect();
    if ids.is_empty() {
        return "original".into();
    }
    let joined = ids.join("+");
    if joined.len() > 100 {
        format!("n{}_{}", ids.len(), ids.last().unwrap())
    } else {
        joined
    }
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
    sets.into_iter()
        .filter(|a| seen.insert(a.clone()))
        .map(|active| OutSet { name: name_of(plan, &active), active })
        .collect()
}
