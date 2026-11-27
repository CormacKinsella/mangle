//! FASTA reading and (multi-target, line-wrapped) writing.

use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

/// Whole reference held in memory (case preserved). Outputs are streamed, so peak memory is
/// roughly one copy of the input plus the largest single variant region.
#[derive(Clone, Debug, Default)]
pub struct Reference {
    pub names: Vec<String>,
    pub seqs: Vec<Vec<u8>>,
}

impl Reference {
    pub fn from_pairs(pairs: &[(&str, &str)]) -> Reference {
        Reference {
            names: pairs.iter().map(|p| p.0.to_string()).collect(),
            seqs: pairs.iter().map(|p| p.1.as_bytes().to_vec()).collect(),
        }
    }
}

/// Longest single line accepted (decimal 2 GB). Unwrapped chromosomes are legitimate, but a file with no
/// line breaks at all (e.g. a binary file renamed .fa) must not be read into memory as one line.
pub const MAX_LINE_BYTES: usize = 2_000_000_000;

/// Most records (contigs) accepted. Each record costs memory beyond its bytes, so a file of millions of tiny
/// records would otherwise use many times its own size.
pub const MAX_RECORDS: usize = 1_000_000;

pub fn read_fasta(path: &Path) -> Result<Reference> {
    read_fasta_with_limits(path, MAX_LINE_BYTES, MAX_RECORDS)
}

pub fn read_fasta_with_limits(path: &Path, max_line: usize, max_records: usize) -> Result<Reference> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut r = BufReader::with_capacity(1 << 20, f);
    let mut refs = Reference::default();
    let mut line = Vec::new();
    let mut line_no = 0usize;
    loop {
        line.clear();
        line_no += 1;
        if r.by_ref().take(max_line as u64 + 1).read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.len() > max_line {
            bail!("{}: a line is longer than {} bytes; is this a FASTA file?", path.display(), max_line);
        }
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if line.first() == Some(&b'>') {
            let name = String::from_utf8_lossy(&line[1..]);
            let name = name.split_whitespace().next().unwrap_or("").to_string();
            if name.is_empty() {
                bail!("FASTA record with empty name");
            }
            // fail at the first empty record instead of after reading the whole file
            if refs.seqs.last().is_some_and(|q| q.is_empty()) {
                bail!("{}: sequence '{}' is empty (a header with no sequence before the next header)", path.display(), refs.names.last().unwrap());
            }
            if refs.names.len() >= max_records {
                bail!("{}: more than {max_records} sequences; is this a FASTA file?", path.display());
            }
            refs.names.push(name);
            refs.seqs.push(Vec::new());
        } else if !line.is_empty() {
            match refs.seqs.last_mut() {
                Some(s) => s.extend(line.iter().filter(|b| !b.is_ascii_whitespace())),
                None => bail!(
                    "{}: line {line_no} is not a valid FASTA header: the file must start with '>' followed by a name, but sequence data came first. Is this a FASTA file?",
                    path.display()
                ),
            }
        }
    }
    if refs.seqs.is_empty() {
        bail!("no sequences in {}", path.display());
    }
    if let Some(i) = refs.seqs.iter().position(|s| s.is_empty()) {
        bail!("sequence '{}' is empty", refs.names[i]);
    }
    Ok(refs)
}

/// Writes wrapped FASTA to several sinks at once (per-sample file + optional combined file).
pub struct FastaWriter<'a> {
    outs: Vec<Box<dyn Write + 'a>>,
    col: usize,
    width: usize,
}

impl<'a> FastaWriter<'a> {
    pub fn new(outs: Vec<Box<dyn Write + 'a>>, width: usize) -> Self {
        FastaWriter { outs, col: 0, width: width.max(1) }
    }

    pub fn begin(&mut self, name: &str) -> Result<()> {
        for o in &mut self.outs {
            writeln!(o, ">{name}")?;
        }
        self.col = 0;
        Ok(())
    }

    pub fn seq(&mut self, mut s: &[u8]) -> Result<()> {
        while !s.is_empty() {
            let n = (self.width - self.col).min(s.len());
            for o in &mut self.outs {
                o.write_all(&s[..n])?;
                if self.col + n == self.width {
                    o.write_all(b"\n")?;
                }
            }
            self.col = (self.col + n) % self.width;
            s = &s[n..];
        }
        Ok(())
    }

    pub fn end(&mut self) -> Result<()> {
        if self.col > 0 {
            for o in &mut self.outs {
                o.write_all(b"\n")?;
            }
        }
        self.col = 0;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        for o in &mut self.outs {
            o.flush()?;
        }
        Ok(())
    }
}
