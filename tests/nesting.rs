use mangle::apply::{complement, render_contig_vec, revcomp};
use mangle::io::Reference;
use mangle::plan::resolve;
use mangle::sets::{Mode, build_sets};
use mangle::spec::parse_spec;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

fn random_ref(len: usize, seed: u64) -> Reference {
    let mut r = ChaCha8Rng::seed_from_u64(seed);
    let s: String = (0..len).map(|_| b"ACGT"[r.gen_range(0..4)] as char).collect();
    Reference::from_pairs(&[("chr1", &s)])
}

const INV_SNP: &str = "id type length parent\nvar1 inversion 40 .\nvar2 snp . var1\n";

#[test]
fn user_example_gives_three_then_two() {
    let refs = random_ref(2000, 1);
    let plan = resolve(&parse_spec(INV_SNP).unwrap(), &refs, 7, 10).unwrap();
    let all = build_sets(&plan, Mode::Cumulative, false);
    let names: Vec<_> = all.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(names, ["product_0000", "product_0001", "product_0002"]);
    let dropped = build_sets(&plan, Mode::Cumulative, true);
    assert_eq!(dropped.len(), 2);
    assert_eq!(dropped[0].name, "product_0001");
}

#[test]
fn original_is_identical_to_input() {
    let refs = random_ref(2000, 2);
    let plan = resolve(&parse_spec(INV_SNP).unwrap(), &refs, 3, 10).unwrap();
    let original = &build_sets(&plan, Mode::Cumulative, false)[0];
    assert_eq!(render_contig_vec(&plan, &refs, &original.active, 0), refs.seqs[0]);
}

#[test]
fn snp_inside_inversion_is_complemented_at_mirrored_position() {
    for seed in 0..50 {
        let refs = random_ref(1500, 100 + seed);
        let plan = resolve(&parse_spec(INV_SNP).unwrap(), &refs, seed, 10).unwrap();
        let sets = build_sets(&plan, Mode::Cumulative, false);
        let inv_only = render_contig_vec(&plan, &refs, &sets[1].active, 0);
        let both = render_contig_vec(&plan, &refs, &sets[2].active, 0);
        let (inv, snp) = (&plan.nodes[0], &plan.nodes[1]);
        let seq = &refs.seqs[0];

        // inversion alone
        let mut expect = seq.clone();
        expect.splice(inv.off..inv.off + inv.len, revcomp(&seq[inv.off..inv.off + inv.len]));
        assert_eq!(inv_only, expect, "seed {seed}");

        // the SNP is inside the inversion, in the parent's frame
        assert!(snp.off < inv.len);
        assert_eq!(snp.ref_pos, inv.off + snp.off);
        assert_eq!(snp.refb, seq[snp.ref_pos]);
        assert_ne!(snp.alt, snp.refb);

        // outputs differ at exactly one position: the mirrored one, holding complement(alt)
        assert_eq!(both.len(), inv_only.len());
        let diffs: Vec<usize> = (0..both.len()).filter(|&i| both[i] != inv_only[i]).collect();
        let mirrored = inv.off + inv.len - 1 - snp.off;
        assert_eq!(diffs, vec![mirrored], "seed {seed}");
        assert_eq!(both[mirrored], complement(snp.alt));
    }
}

#[test]
fn same_seed_same_plan_different_seed_differs() {
    let refs = random_ref(5000, 9);
    let spec = parse_spec("a inv 100-300 .\nb snp . a\nc del 50 .\nd ins 20 .\ne cnv 30 .\n").unwrap();
    let p1 = resolve(&spec, &refs, 5, 10).unwrap();
    let p2 = resolve(&spec, &refs, 5, 10).unwrap();
    let p3 = resolve(&spec, &refs, 6, 10).unwrap();
    let key = |p: &mangle::plan::Plan| p.nodes.iter().map(|n| (n.off, n.len, n.alt)).collect::<Vec<_>>();
    assert_eq!(key(&p1), key(&p2));
    assert_ne!(key(&p1), key(&p3));
}

#[test]
fn variant_position_is_stable_across_outputs() {
    // the whole point of resolving the plan once: var1 must be at the same place in every output
    let refs = random_ref(4000, 4);
    let spec = parse_spec("a del 30 .\nb inv 60 .\nc snp . .\n").unwrap();
    let plan = resolve(&spec, &refs, 11, 10).unwrap();
    let sets = build_sets(&plan, Mode::Cumulative, false);
    let a = &plan.nodes[0];
    let one = render_contig_vec(&plan, &refs, &sets[1].active, 0);
    let all = render_contig_vec(&plan, &refs, &sets[3].active, 0);
    // everything before the first (lowest-coordinate) variant is untouched in both
    let first = plan.nodes.iter().map(|n| n.off).min().unwrap();
    assert_eq!(one[..first], refs.seqs[0][..first]);
    assert_eq!(all[..first], refs.seqs[0][..first]);
    assert_eq!(one.len(), refs.seqs[0].len() - a.len);
}

#[test]
fn nested_children_always_fit_inside_parent() {
    for seed in 0..100 {
        let refs = random_ref(3000, seed);
        let spec = parse_spec("p inv 30 .\nc1 snp . p\nc2 del 10 p\nc3 ins 5 p\nc4 snp . p\n").unwrap();
        let plan = resolve(&spec, &refs, seed, 10).unwrap();
        let p = &plan.nodes[0];
        let mut spans: Vec<(usize, usize)> = (1..5)
            .map(|i| (plan.nodes[i].off, plan.nodes[i].off + plan.nodes[i].span()))
            .collect();
        for &(s, e) in &spans {
            assert!(e <= p.len && s <= e);
        }
        // no overlap among sibling consuming spans
        spans.sort();
        for w in spans.windows(2) {
            assert!(w[0].1 <= w[1].0, "seed {seed}: overlap {:?}", w);
        }
    }
}

#[test]
fn children_too_big_for_parent_is_an_error() {
    let refs = random_ref(1000, 1);
    let spec = parse_spec("p inv 10 .\nc del 50 p\n").unwrap();
    assert!(resolve(&spec, &refs, 1, 10).is_err());
}

#[test]
fn cnv_with_nested_snp_modifies_one_copy() {
    let refs = random_ref(2000, 6);
    let spec = parse_spec("u cnv 25 . 4\ns snp . u\n").unwrap();
    let plan = resolve(&spec, &refs, 2, 10).unwrap();
    let sets = build_sets(&plan, Mode::Cumulative, false);
    let u = &plan.nodes[0];
    let out = render_contig_vec(&plan, &refs, &sets[2].active, 0);
    assert_eq!(out.len(), refs.seqs[0].len() + 3 * u.len);
    let unit = &refs.seqs[0][u.off..u.off + u.len];
    for c in 0..3 {
        assert_eq!(&out[u.off + c * u.len..u.off + (c + 1) * u.len], unit);
    }
    let last = &out[u.off + 3 * u.len..u.off + 4 * u.len];
    assert_eq!(last.iter().zip(unit).filter(|(a, b)| a != b).count(), 1);
}

#[test]
fn snp_inside_insertion_lands_in_inserted_sequence() {
    let refs = random_ref(2000, 8);
    let spec = parse_spec("i ins 30 .\ns snp . i\n").unwrap();
    let plan = resolve(&spec, &refs, 3, 10).unwrap();
    let sets = build_sets(&plan, Mode::Cumulative, false);
    let i = &plan.nodes[0];
    let s = &plan.nodes[1];
    let ins_only = render_contig_vec(&plan, &refs, &sets[1].active, 0);
    let both = render_contig_vec(&plan, &refs, &sets[2].active, 0);
    assert_eq!(ins_only.len(), refs.seqs[0].len() + 30);
    assert_eq!(&ins_only[i.off..i.off + 30], &i.gen_seq[..]);
    let g = s.gen_off.unwrap();
    assert_eq!(both[i.off + g], s.alt);
    assert_eq!(s.refb, i.gen_seq[g]);
}

#[test]
fn case_is_preserved_through_inversion() {
    let refs = Reference::from_pairs(&[("c", &"acgtNNACGTacgt".repeat(20))]);
    let spec = parse_spec("v inv 10 .\n").unwrap();
    let plan = resolve(&spec, &refs, 1, 0).unwrap();
    let sets = build_sets(&plan, Mode::Cumulative, false);
    let out = render_contig_vec(&plan, &refs, &sets[1].active, 0);
    assert_eq!(out.len(), refs.seqs[0].len());
    assert!(out.iter().all(|b| b"acgtACGTNn".contains(b)));
}

#[test]
fn single_mode_includes_ancestors() {
    let refs = random_ref(2000, 1);
    let plan = resolve(&parse_spec("a inv 30 .\nb snp . a\nc del 10 .\n").unwrap(), &refs, 1, 10).unwrap();
    let names: Vec<_> = build_sets(&plan, Mode::Single, false).into_iter().map(|s| s.name).collect();
    assert_eq!(names, ["product_0000", "product_0001", "product_0002", "product_0003"]);
}

#[test]
fn variants_avoid_n_runs() {
    let refs = Reference::from_pairs(&[("c", &format!("{}{}{}", "ACGT".repeat(50), "N".repeat(200), "ACGT".repeat(50)))]);
    for seed in 0..50 {
        let plan = resolve(&parse_spec("v inv 20 .\nw del 5 .\n").unwrap(), &refs, seed, 10).unwrap();
        for n in &plan.nodes {
            assert!(!refs.seqs[0][n.off..n.off + n.len.max(1)].contains(&b'N'));
        }
    }
}

#[test]
fn vcf_is_sorted_with_children_at_parent_start() {
    for seed in 0..200 {
        let refs = random_ref(1500, seed);
        let spec = parse_spec("p inv 30 .\nc del 5 p\nd ins 4 p\ne snp . p\nf dup 6 p\n").unwrap();
        let plan = resolve(&spec, &refs, seed, 10).unwrap();
        let all = vec![true; plan.nodes.len()];
        let mut buf = Vec::new();
        mangle::truth::write_vcf(&mut buf, &plan, &refs, &all, "t").unwrap();
        let pos: Vec<usize> = String::from_utf8(buf)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
            .collect();
        assert!(pos.windows(2).all(|w| w[0] <= w[1]), "seed {seed}: {pos:?}");
        assert!(pos.iter().all(|&p| p >= 1));
    }
}

#[test]
fn vcf_omits_variants_inside_insertions() {
    let refs = random_ref(1500, 3);
    let plan = resolve(&parse_spec("i ins 30 .\ns snp . i\n").unwrap(), &refs, 1, 10).unwrap();
    let mut buf = Vec::new();
    mangle::truth::write_vcf(&mut buf, &plan, &refs, &[true, true], "t").unwrap();
    let t = String::from_utf8(buf).unwrap();
    assert!(t.contains("\ti\t") && !t.contains("\ts\t"));
}

#[test]
fn bed_matches_plan_and_skips_inserted_children() {
    let refs = random_ref(2000, 5);
    let plan = resolve(&parse_spec("p inv 30 .\nc snp . p\ni ins 10 .\nj snp . i\n").unwrap(), &refs, 2, 10).unwrap();
    let mut buf = Vec::new();
    mangle::truth::write_bed(&mut buf, &plan, &refs).unwrap();
    let t = String::from_utf8(buf).unwrap();
    assert_eq!(t.lines().count(), 3, "{t}");
    assert!(t.contains(&format!("chr1\t{}\t{}\tp\t", plan.nodes[0].off, plan.nodes[0].off + 30)));
    let i = &plan.nodes[2];
    assert!(t.contains(&format!("chr1\t{}\t{}\ti\t", i.off, i.off))); // zero-length insertion
    assert!(!t.contains("\tj\t"));
}
