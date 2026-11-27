use mangle::apply::{contig_len, render_contig_vec};
use mangle::io::{FastaWriter, Reference};
use mangle::limits::{OutputShape, check_output_size, estimate_output_bytes};
use mangle::plan::resolve;
use mangle::sets::{Mode, build_sets};
use mangle::spec::{parse_size, parse_spec};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::process::Command;

fn random_ref(len: usize, seed: u64) -> Reference {
    let mut r = ChaCha8Rng::seed_from_u64(seed);
    let s: String = (0..len).map(|_| b"ACGT"[r.gen_range(0..4)] as char).collect();
    Reference::from_pairs(&[("chr1", &s), ("chr2", &s[..len / 2])])
}

// every type, with nesting inside inv, dup, cnv and ins
const SPEC: &str = "id type length parent copies\n\
    a inv 60 .\n a1 snp . a\n a2 ins 7 a\n a3 dup 5 a\n\
    b cnv 30 . 4\n b1 snp . b\n b2 del 4 b\n\
    c ins 25 .\n c1 inv 10 c\n c2 snp . c\n\
    d dup 40 .\n e del 20 .\n";

fn shape(h: usize, combined: bool) -> OutputShape {
    OutputShape { haplotypes_per_product: h, combined, line_width: 80 }
}

#[test]
fn computed_lengths_match_rendered_sequence() {
    for seed in 0..30 {
        let refs = random_ref(3000, seed);
        let plan = resolve(&parse_spec(SPEC).unwrap(), &refs, seed, 10).unwrap();
        for mode in [Mode::Cumulative, Mode::Single] {
            for s in build_sets(&plan, mode, false) {
                for c in 0..2 {
                    assert_eq!(contig_len(&plan, &refs, &s.active, c), render_contig_vec(&plan, &refs, &s.active, c).len(), "seed {seed} {}", s.name);
                }
            }
        }
    }
}

#[test]
fn estimate_equals_bytes_actually_written() {
    for width in [1usize, 7, 80, 100_000] {
        let refs = random_ref(3000, 5);
        let plan = resolve(&parse_spec(SPEC).unwrap(), &refs, 5, 10).unwrap();
        let sets = build_sets(&plan, Mode::Cumulative, false);
        let mut written = 0u64;
        for s in &sets {
            let mut buf: Vec<u8> = Vec::new();
            {
                let mut w = FastaWriter::new(vec![Box::new(&mut buf)], width);
                for (ci, name) in ["chr1", "chr2"].iter().enumerate() {
                    w.begin(&format!("{}#0#{}", s.name, name)).unwrap();
                    w.seq(&render_contig_vec(&plan, &refs, &s.active, ci)).unwrap();
                    w.end().unwrap();
                }
                w.flush().unwrap();
            }
            written += buf.len() as u64;
        }
        let sh = OutputShape { haplotypes_per_product: 1, combined: false, line_width: width };
        assert_eq!(estimate_output_bytes(&plan, &refs, &sets, &sh), written, "width {width}");
    }
}

#[test]
fn drop_input_haplotypes_and_combined_change_the_estimate() {
    let refs = random_ref(3000, 9);
    let plan = resolve(&parse_spec(SPEC).unwrap(), &refs, 9, 10).unwrap();
    let with = build_sets(&plan, Mode::Single, false);
    let without = build_sets(&plan, Mode::Single, true);
    let (e_with, e_without) = (estimate_output_bytes(&plan, &refs, &with, &shape(1, false)), estimate_output_bytes(&plan, &refs, &without, &shape(1, false)));
    assert!(e_with > e_without, "the unchanged input is counted unless dropped");
    assert_eq!(estimate_output_bytes(&plan, &refs, &with, &shape(3, false)), 3 * e_with);
    assert_eq!(estimate_output_bytes(&plan, &refs, &with, &shape(1, true)), 2 * e_with);
}

#[test]
fn limit_is_enforced() {
    let refs = random_ref(3000, 9);
    let plan = resolve(&parse_spec(SPEC).unwrap(), &refs, 9, 10).unwrap();
    let sets = build_sets(&plan, Mode::Cumulative, false);
    let est = estimate_output_bytes(&plan, &refs, &sets, &shape(1, false));
    assert_eq!(check_output_size(&plan, &refs, &sets, &shape(1, false), est).unwrap(), est);
    let err = check_output_size(&plan, &refs, &sets, &shape(1, false), est - 1).unwrap_err().to_string();
    assert!(err.contains("--max_output_size"), "{err}");
    assert!(check_output_size(&plan, &refs, &sets, &shape(2, false), est).is_err());
}

#[test]
fn size_suffixes_and_overflow() {
    assert_eq!(parse_size("10g").unwrap(), 10_000_000_000);
    assert_eq!(parse_size("500M").unwrap(), 500_000_000);
    assert_eq!(parse_size("2k").unwrap(), 2000);
    assert_eq!(parse_size("1_000").unwrap(), 1000);
    assert!(parse_size("99999999999999999999g").is_err());
    assert!(parse_size("20000000000g").is_err());
    assert!(parse_size("abc").is_err());
}

#[test]
fn cli_refuses_before_writing_anything() {
    let dir = std::env::temp_dir().join(format!("mangle_limits_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fa = dir.join("in.fa");
    let spec = dir.join("spec.tsv");
    std::fs::write(&fa, format!(">chr1\n{}\n", "ACGT".repeat(2500))).unwrap();
    std::fs::write(&spec, "v1 inv 100 .\nv2 del 50 .\n").unwrap();
    let out = dir.join("out");
    let run = |limit: &str| {
        Command::new(env!("CARGO_BIN_EXE_mangle"))
            .args(["-s", spec.to_str().unwrap(), "-o", out.to_str().unwrap(), "--max_output_size", limit, fa.to_str().unwrap()])
            .output()
            .unwrap()
    };
    let bad = run("10k"); // 3 sequences of ~10 kb
    assert!(!bad.status.success());
    assert!(String::from_utf8_lossy(&bad.stderr).contains("max_output_size"));
    assert!(!out.exists(), "nothing may be written when refused");
    assert!(run("1m").status.success());
    assert!(out.join("original.fa").exists());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn oversized_variants_are_rejected_regardless_of_seed() {
    // growth = len for ins/dup, len x (copies-1) for cnv; limit 1 GB, cnv copies <= 500
    assert!(parse_spec("a cnv 100m . 11\n").is_ok()); // exactly 1 GB
    assert!(parse_spec("a cnv 100m . 12\n").is_err());
    assert!(parse_spec("a cnv 10 . 501\n").is_err());
    assert!(parse_spec("a ins 2g .\n").is_err());
    assert!(parse_spec("a dup 2g .\n").is_err());
    assert!(parse_spec("a ins 1k-2g .\n").is_err(), "the top of a range counts, so the seed can't matter");
    assert!(parse_spec("a inv 5g .\n").is_ok(), "types that don't grow are not limited here");
}
