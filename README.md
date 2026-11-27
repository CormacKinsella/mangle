[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)
[![Rust 1.100](https://img.shields.io/badge/rust-1.100-orange?logo=rust)](https://www.rust-lang.org)

<!-- logo -->
<img src="assets/logo.svg" alt="mangle dotplot logo" width="480">

## At a glance

`mangle` is a support tool for pangenome applications that simulates structural variation in DNA sequences.

It simulates `SNPs`, `indels`, `inversions`, `tandem duplications`, and `CNVs` as either `top-level` or `nested` variants (i.e., a variant within a variant (within a variant (within a variant))).

Unlike many simulators that produce a single mutated sequence, `mangle` outputs a **family of related sequences** that are ready for graph construction.

Note that each product is based on the full input, so output can be over `n + 1` sequences in size. Keep this in mind for large input files and/or many variants. mangle estimates the total output size before writing anything and refuses to run above `--max_output_size` (default `10g`).

## Quick start

The Rust toolchain lives in the repo software environment, so [pixi](https://pixi.prefix.dev/latest/) is the only prerequisite. 

Paths are relative to where you run it, so use the repo root or absolute paths:

```
pixi run mangle --help
pixi run mangle input.fa --spec spec.tsv
```

`pixi run mangle` depends on the `build` task (i.e. `cargo build --release`).

`-s` or `--spec` is used to provide a variant catalogue. The [file](#variant-catalogue) should be a four-column `TSV` (the fifth column is optional):

```
id    type       length   parent    copies
seq1  inversion  1000     .
seq2  snp        1        seq1
seq3  cnv        15       .         5

# Note: in the parent column a full stop character denotes a top-level variant (no nesting)
```

### Full installation

```
pixi run cargo install --path . --root ~/.local
export PATH="$HOME/.local/bin:$PATH"
mangle --version
```

Re-run the `cargo install` line after changing the source to update the installed copy. The built binary does not need pixi, so it can be copied to a cluster as a single file (same Linux architecture and a compatible glibc).

## Quick visualisation of outputs

The tool [dnadotplot](https://github.com/quadram-institute-bioscience/dnadotplot) is provided in the `pixi` environment along with tasks to run it.

```
# all pairwise plots of the products in outputs/ -> plots/<A>_vs_<B>.png
pixi run plots [input_dir=outputs] [output_dir=plots] [window_size=7]

# one pair; only the first record of each file is plotted
pixi run plotspecific <a.fa> <b.fa> <out.png> [window]
```

## Variant catalogue

A `TSV` file describing the variants to simulate, one variant per line. The header line is optional, and `# comment lines` are ignored.

Column headers: `id type length parent [copies]`

| column | meaning |
|---|---|
| `id` | Unique variant/product sequence identifier (`[A-Za-z0-9_.-]`, cannot be `original`) |
| `type` | Variant type (`snp` `ins` `del` `inv` `dup` `cnv`, long names like `inversion` work, and `snv` for `snp`) |
| `length` | Variant length (exact or a range: `500`, `200-400`/`200..400`, suffixes `k/m/g`) |
| `parent` | id of an enclosing variant, or `.`, `-`, `0`, `none` for top-level. Parents must be declared first; `snp` and `del` cannot be parents |
| `copies` | total copy number of the region, applies to `cnv` only (default 3, minimum 2, maximum 500), `dup` is always 2 |

What each type does: `snp` substitutes one base; `ins` inserts a random `ACGT` sequence of the given length; `del` removes
the region; `inv` reverse-complements it; `dup` makes one extra tandem copy (2 copies in total); `cnv` makes `copies` tandem
copies in total.

## Run mode options

`--mode` decides which combinations of the catalogue's variants become output sequences. Every output sequence contains a subset of a resolved variant catalogue (all random choices are drawn once), so a given variant will sit at the same place in every sequence that contains it, whichever mode you pick.

### `cumulative` (default)

Cumulative addition of variants, in catalogue order. It outputs the unmodified input first, then the first catalogue line, e.g., `var1`, then `var1+var2`, then `var1+var2+var3`, etc.

For a catalogue of `n` variants you get `n` products (plus `original`).

### `single`

Each variant on its own, in catalogue order. Since a nested variant needs its parent, these outputs will contain all the respective ancestors (a child of a child includes both the parent and the grandparent). Siblings and unrelated top-level variants are not combined.

For a catalogue of `n` variants you get `n` products (plus `original`).

#### Example

```
id    type       length  parent
var1  inversion  1000    .
var2  snp        .       var1
var3  del        200     .
```

| `--mode` | sequences written |
|---|---|
| `cumulative` | `original`, `var1`, `var1+var2`, `var1+var2+var3` |
| `single` | `original`, `var1`, `var1+var2`, `var3` |

## Output size limit

Before writing anything, mangle computes the exact size of the FASTA it is about to produce: every sequence in the chosen `--mode`, including the unchanged input unless `--drop_input` is set, and doubled when `--combined` also writes `all.fa`. If that exceeds `--max_output_size` (default `10g`, decimal `k/m/g` suffixes) it stops with an error and writes nothing. Raise the limit explicitly (e.g. `--max_output_size 200g`) when you really mean it. Truth files (BED, VCF, catalogue) are small and not counted. The estimate is printed when the run starts and saved in `mangle.log`.

Separately, single variants are limited in how much they can add to the sequence, because each one is built in memory: `ins` and `dup` add their length, `cnv` adds length × (copies − 1), and the limit is 1 GB each (for a length range the top of the range is checked). These two limits (500 copies, 1 GB) are fixed. A spec may also have at most 100,000 variants, and before placing anything mangle checks that the top-level variants (plus `--min_gap`) could fit in the non-N bases of the reference at all.

## Outputs

| file | when | contents |
|---|---|---|
| `NAME.fa` | always | one per product; the unmodified input is `original.fa` |
| `resolved_catalogue.tsv` | always | your variant catalogue after the RNG has acted: one row per variant with its contig, position, length, parent, and for SNPs the ref/alt bases. This is the only file listing variants inside inserted sequence (`ins_offset`) |
| `mangle.log` | always | exact command line, parameters, warnings, md5 of the variant catalogue |
| `variants.bed` | always | every variant that has an original-reference coordinate, 0-based; insertions are zero-length. Variants inside inserted sequence are not listed |
| `all.fa` | only with `--combined` | every output sequence (including `original`, unless `--drop_input`) in one multi-sample FASTA |
| `NAME.vcf` | only with `--vcf` | truth VCF for each product |

See "Sequence names" below for the FASTA record naming.

## Sequence names

Every FASTA record is named `PRODUCT#0#CONTIG` ([PanSN](https://github.com/pangenome/PanSN-spec)) for compatibility with graph building tools.

* **PRODUCT** is `original` for the unmodified input, otherwise the ids of the variants it contains, in
  catalogue order, joined with `+`.
* **0** is the haplotype field. It is always `0` for now.
* **CONTIG** is the input name up to the first space.

| input header | product | output header |
|---|---|---|
| `>chr1 Homo sapiens chromosome 1` | the input | `>original#0#chr1` |
| `>chr1 Homo sapiens chromosome 1` | `var1` | `>var1#0#chr1` |
| `>scaffold_2` | `var1` + nested `var2` | `>var1+var2#0#scaffold_2` |
| `>HG002#1#chr1` (already PanSN) | `var1` | `>var1#0#chr1` (only the last `#` field is kept) |

* `--mode single` names a nested variant with its ancestors, since it needs them: `var2` inside `var1` is
  `var1+var2`, and an independent `var3` is just `var3`.
* If the `+`-joined name would exceed 100 characters it is shortened to `n<COUNT>_<LAST_ID>`, for example
  `n7_a_rather_long_variant_identifier_07` for a product with seven variants. In the default cumulative mode
  that means the first 7 variants of the catalogue, ending in the named one. The shortened name is also what
  `mangle.log` lists, so the mapping is not recorded anywhere else yet.
* Input contigs whose names collide after the PanSN prefix is stripped (`HG002#1#chr1` and `HG002#2#chr1`)
  are an error, since the output names would not be unique. Rename them first.

## Column reference

**`variants.bed`** has no header line (standard BED, so bedtools/IGV read it directly). Tab-separated,
BED6 plus two extra columns, 0-based half-open, in the coordinates of the **original input**:

| # | column | meaning |
|---|---|---|
| 1 | chrom | input contig name (up to the first space) |
| 2 | start | 0-based start |
| 3 | end | end, exclusive. Insertions are zero-length (`start == end`) |
| 4 | name | variant `id` from the catalogue |
| 5 | score | always `0` (placeholder required by BED6) |
| 6 | strand | always `.` |
| 7 | type | `snp` `ins` `del` `inv` `dup` `cnv` |
| 8 | parent | `id` of the enclosing variant, or `.` if top-level |

A nested variant is placed in its parent's unmodified frame, so it falls inside the parent's interval. One
row per variant, whether or not it appears in a given product, except variants inside inserted sequence: they have no
original coordinate and are not listed here (see `resolved_catalogue.tsv`).

**`resolved_catalogue.tsv`** has a header line:

| column | meaning |
|---|---|
| `id`, `type`, `parent` | as in the variant catalogue (`parent` is `.` if top-level) |
| `contig` | input contig name |
| `ref_pos_1based` | 1-based start in the original input. For a variant inside inserted sequence, the insertion point |
| `len` | drawn length (1 for `snp`; inserted length for `ins`; unit length for `dup`/`cnv`) |
| `copies` | total copy number (`dup` 2, `cnv` n, otherwise 1) |
| `ref`, `alt` | `snp` only: reference and substituted base (`.` otherwise) |
| `ins_offset` | for a variant inside inserted sequence, its 0-based offset within that sequence, else `.` |

For `ins`/`del`, `ref_pos_1based` is one more than the VCF POS, because the VCF uses the anchor-base
convention.

## Nesting semantics

* A child is placed inside its parent's input region; siblings never overlap. Positions are
  independent of which variants end up active.
* Child of `inv`: applied before the inversion, so it appears reverse-complemented
  (SNP at the mirrored position, complemented allele).
* Child of `dup`/`cnv`: lands in one extra copy (the last); the other copies stay pristine.
* Child of `ins`: placed inside the generated inserted sequence.
* `snp` and `del` cannot be parents.
* Top-level variants avoid runs of `N`. Case (soft-masking) is preserved.

## Truth VCF

Coordinates are in the **original reference frame**. A nested child is reported in its parent's unmodified frame
(a SNP inside an inversion is reported where it sits on the reference, not where it ends up after
reverse-complementing). `ins`/`del` use the VCF anchor-base convention, so their POS is one less than
`ref_pos_1based` in `resolved_catalogue.tsv`. Variants *inside inserted sequence* have no reference coordinate and are
left out of the VCF (they are in `resolved_catalogue.tsv`, with `ins_offset`). The VCF is sites-only (no sample
columns). SNPs carry real REF/ALT bases; the other types use symbolic ALTs (`<DEL>`, `<INS>`, `<INV>`, `<DUP>`,
`<CNV>`) with `SVTYPE`, `SVLEN`, `END`, `CN` (dup/cnv) and `PARENT` (nested) in INFO. A spot check with
`bcftools norm -f ref.fa -c e` (bcftools is in the pixi environment) reported no errors, but the test suite does not run it.
Not yet provided: coordinates in each output sequence's own frame.

## Options reference

```
Options:
  -s, --spec <SPEC>
          Variant catalogue (TSV: id type length parent [copies])

  -o, --outdir <OUTDIR>
          Output directory
          
          [default: outputs]

      --seed <SEED>
          RNG seed (same seed + same inputs => identical output)
          
          [default: 42]

      --mode <MODE>
          How variants are combined into output sequences

          Possible values:
          - cumulative: Ladder: original, v1, v1+v2, v1+v2+v3, ... (spec order; parents always precede children)
          - single:     Each variant on its own (plus its ancestors, since a nested variant needs its parent)
          
          [default: cumulative]

      --drop_input
          Do not emit the input sequence, only the products

      --vcf
          Also write a truth VCF for each product

      --combined
          Also write every output into one multi-sample FASTA (all.fa), e.g. for PGGB

      --min_gap <MIN_GAP>
          Minimum distance between top-level variants (bp)
          
          [default: 10]

      --max_output_size <MAX_OUTPUT_SIZE>
          Refuse to run if the FASTA output would exceed this size (suffixes k/m/g, decimal)
          
          [default: 10g]

      --line_width <LINE_WIDTH>
          FASTA line width
          
          [default: 80]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

## AI disclosure

`mangle` was written using Claude Code with iterative feedback. 
