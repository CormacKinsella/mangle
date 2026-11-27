# Changelog

## [0.2.0] - unreleased

### Added
- `--max_output_size` (default `10g`, suffixes `k/m/g`): mangle works out how large the total FASTA output will be before writing anything, and stops with a clear message if it is over the limit. It counts every sequence written, including the unchanged input (unless `--drop_input`) and the extra `all.fa` copy from `--combined`. The estimate is printed at run time and recorded in `mangle.log`.
- Per-variant size limits: a `cnv` may have at most 500 copies, and no single variant (`ins`, `dup`, `cnv`) may add more than 1 GB to the sequence.
- A test that pins the exact output for a fixed seed, so a dependency upgrade that changes what a seed produces is caught.
- Tests for the size estimate (it matches the bytes actually written) and for the limit itself.

### Fixed
- Very large sizes in the catalogue (for example `20000000000g`) now give an error instead of silently wrapping around.

## [0.1.0] - 2026-11-27

**Full Changelog**: https://github.com/CormacKinsella/mangle/commits/0.1.0

# mangle v0.1.0

First release. `mangle` simulates structural variation and writes a family of related haplotypes for pangenome construction and testing.

## What's in it
- Variant types: `snp`, `ins`, `del`, `inv`, `dup`, `cnv`, as top-level or nested variants.
- Inputs: a FASTA (`.fa`/`.fasta`/`.fna`) and a variant catalogue TSV (`id type length parent [copies]`). Lengths can be fixed or a range.
- Modes: `cumulative` (default, `var1`, `var1+var2`, ...) and `single` (each variant with its ancestors).
- Outputs: one FASTA per product (PanSN names `product#0#contig`), `variants.bed`, `resolved_catalogue.tsv`, `mangle.log`, and optional `--vcf` and `--combined`.
- Reproducible: the same seed and inputs give identical output. N stretches are never touched. Software environment provided via lock files

## How it works (input to output)
1. **Parse** (`spec.rs`): the TSV becomes a validated catalogue. Parents must be declared first, and `snp` and `del` cannot be parents.
2. **Plan** (`plan.rs`): the seeded RNG draws every position, length, allele and inserted sequence once. A top-level variant lives in reference coordinates. A nested variant lives at an offset inside its parent's input region. Top-levels may touch but never overlap. N regions are avoided.
3. **Select** (`sets.rs`): the mode picks which subsets of the plan become products. Every product is a subset of the same plan, so a variant sits at the same place in every product that contains it.
4. **Apply** (`apply.rs`): each product is built recursively. A node (variant) takes its input region, swaps in each active child's realised output, then applies its own operation. A SNP inside an inversion is therefore reverse-complemented with it.
5. **Write** (`io.rs`, `truth.rs`): the reference is held in memory and outputs are streamed to FASTA. The BED, resolved catalogue, log and optional VCFs are written alongside. Truth coordinates are in the original reference frame.

Memory is about 1x the input. Disk is N+1 genome copies for N variants in current default modes.
