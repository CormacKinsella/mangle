//! Pins the exact output for one seed. If this fails after a dependency bump (e.g. `rand`), the same seed
//! no longer gives the same results: that breaks the reproducibility promise, so treat it as a breaking change.
use mangle::io::Reference;
use mangle::plan::resolve;
use mangle::spec::parse_spec;

#[test]
fn seed_42_gives_pinned_plan() {
    // reference built without the RNG, so only the plan's randomness is under test
    let seq: String = (0..3000usize).map(|i| b"ACGT"[(i * i / 7 + i) % 4] as char).collect();
    let refs = Reference::from_pairs(&[("chr1", &seq)]);
    let spec = parse_spec("a inv 100-300 .\nb snp . a\nc del 50 .\nd ins 20 .\n").unwrap();
    let plan = resolve(&spec, &refs, 42, 10).unwrap();
    let got: Vec<String> = plan
        .nodes
        .iter()
        .map(|n| format!("{} off={} len={} alt={} ins={}", n.id, n.off, n.len, n.alt as char, String::from_utf8_lossy(&n.gen_seq)))
        .collect();
    assert_eq!(
        got,
        [
            "a off=2590 len=237 alt=N ins=",
            "b off=54 len=1 alt=G ins=",
            "c off=1337 len=50 alt=N ins=",
            "d off=2893 len=20 alt=N ins=CCTTGGGACAACTTATATTA",
        ]
    );
}
