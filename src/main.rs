use anyhow::{Context, Result, bail};
use clap::Parser;
use mangle::io::{FastaWriter, read_fasta};
use mangle::sets::{Mode, build_sets};
use mangle::{apply, plan, spec, truth};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// Simulate nested variation (SNPs, indels, inversions, dups, CNVs) as a family of related sequences.
///
/// Given a variant catalogue, positions are drawn by a seeded RNG; as are lengths when provided as a range (e.g. 200-400), e.g.:
///
///   id    type       length  parent
///   var1  inv        1000    .
///   var2  snp        1       var1
///
/// gives three sequences by default: original, var1, var1+var2.
#[derive(Parser, Debug)]
#[command(version, about, verbatim_doc_comment, before_help = concat!(env!("CARGO_PKG_NAME"), " ", env!("CARGO_PKG_VERSION")))]
struct Args {
    /// Input FASTA (one or more records)
    input: PathBuf,

    /// Variant catalogue (TSV: id type length parent [copies])
    #[arg(short, long)]
    spec: PathBuf,

    /// Output directory
    #[arg(short, long, default_value = "outputs")]
    outdir: PathBuf,

    /// RNG seed (same seed + same inputs => identical output)
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// How variants are combined into output sequences
    #[arg(long, value_enum, default_value_t = Mode::Cumulative)]
    mode: Mode,

    /// Do not emit the input sequence, only the products
    #[arg(long = "drop_input", aliases = ["drop-input", "drop_original", "drop-original"])]
    drop_input: bool,

    /// Also write a truth VCF for each product
    #[arg(long)]
    vcf: bool,

    /// Also write every output into one multi-sample FASTA (all.fa), e.g. for PGGB
    #[arg(long)]
    combined: bool,

    /// Minimum distance between top-level variants (bp)
    #[arg(long = "min_gap", alias = "min-gap", default_value_t = 10)]
    min_gap: usize,

    /// FASTA line width
    #[arg(long = "line_width", alias = "line-width", default_value_t = 80)]
    line_width: usize,
}

fn main() -> Result<()> {
    let a = Args::parse();
    let text = fs::read_to_string(&a.spec).with_context(|| format!("reading {}", a.spec.display()))?;
    let ext_ok = a
        .input
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| matches!(e.to_ascii_lowercase().as_str(), "fasta" | "fna" | "fa"));
    if !ext_ok {
        bail!("input '{}' must be a FASTA file with extension .fasta, .fna or .fa", a.input.display());
    }
    let (specs, warnings) = spec::parse_spec_warn(&text)?;
    for w in &warnings {
        eprintln!("warning: {w}");
    }
    let refs = read_fasta(&a.input)?;
    // Output records are PanSN `sample#hap#contig`. If the input is already PanSN, keep only its contig field.
    let labels: Vec<String> = refs.names.iter().map(|n| n.rsplit('#').next().unwrap_or(n).to_string()).collect();
    {
        let mut seen = std::collections::HashSet::new();
        for (l, n) in labels.iter().zip(&refs.names) {
            if !seen.insert(l) {
                bail!("input contig names collide once PanSN prefixes are stripped ('{n}' -> '{l}'); rename the input contigs");
            }
        }
    }
    let plan = plan::resolve(&specs, &refs, a.seed, a.min_gap)?;
    let sets = build_sets(&plan, a.mode, a.drop_input);
    if sets.is_empty() {
        bail!("nothing to output");
    }
    fs::create_dir_all(&a.outdir)?;

    // The catalogue with every random choice resolved to concrete coordinates and alleles. Always written.
    truth::write_resolved_catalogue(
        &mut BufWriter::new(File::create(a.outdir.join("resolved_catalogue.tsv"))?),
        &plan,
        &refs,
    )?;

    truth::write_bed(&mut BufWriter::new(File::create(a.outdir.join("variants.bed"))?), &plan, &refs)?;

    // Run log: exact command line, parameters, and md5 of the variant catalogue for provenance.
    {
        let mut l = BufWriter::new(File::create(a.outdir.join("mangle.log"))?);
        let unix = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0);
        writeln!(l, "mangle {} (unix time {unix})", env!("CARGO_PKG_VERSION"))?;
        writeln!(l, "command: {}", std::env::args().collect::<Vec<_>>().join(" "))?;
        writeln!(l, "input: {}", a.input.display())?;
        writeln!(l, "spec: {}", a.spec.display())?;
        writeln!(l, "spec_md5: {:x}", md5::compute(text.as_bytes()))?;
        writeln!(l, "seed: {}", a.seed)?;
        writeln!(l, "mode: {:?}", a.mode)?;
        writeln!(l, "drop_input: {}", a.drop_input)?;
        writeln!(l, "min_gap: {}", a.min_gap)?;
        writeln!(l, "vcf: {}", a.vcf)?;
        writeln!(l, "variants: {}", plan.nodes.len())?;
        writeln!(l, "products: {}", sets.iter().map(|s| s.name.as_str()).collect::<Vec<_>>().join(" "))?;
        for w in &warnings {
            writeln!(l, "warning: {w}")?;
        }
    }

    let mut combined = if a.combined {
        Some(BufWriter::new(File::create(a.outdir.join("all.fa"))?))
    } else {
        None
    };

    eprintln!("seed {} | {} variants | {} output sequences", a.seed, plan.nodes.len(), sets.len());
    for s in &sets {
        let mut outs: Vec<Box<dyn Write + '_>> =
            vec![Box::new(BufWriter::new(File::create(a.outdir.join(format!("{}.fa", s.name)))?))];
        if let Some(c) = combined.as_mut() {
            outs.push(Box::new(c));
        }
        let mut w = FastaWriter::new(outs, a.line_width);
        for (ci, cname) in labels.iter().enumerate() {
            // PanSN naming: sample#haplotype#contig
            w.begin(&format!("{}#0#{}", s.name, cname))?;
            let mut res: Result<()> = Ok(());
            apply::render_contig(&plan, &refs, &s.active, ci, &mut |b| {
                if res.is_ok() {
                    res = w.seq(b);
                }
            });
            res?;
            w.end()?;
        }
        w.flush()?;
        // No VCF for the unmodified input: it has no variants relative to itself.
        if a.vcf && s.active.iter().any(|&x| x) {
            truth::write_vcf(
                &mut BufWriter::new(File::create(a.outdir.join(format!("{}.vcf", s.name)))?),
                &plan,
                &refs,
                &s.active,
                &s.name,
            )?;
        }
        eprintln!("  {}", s.name);
    }
    if let Some(c) = combined.as_mut() {
        c.flush()?;
    }
    Ok(())
}
