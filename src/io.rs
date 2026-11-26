//! FASTA reading and (multi-target, line-wrapped) writing.

use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
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

pub fn read_fasta(path: &Path) -> Result<Reference> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut r = BufReader::with_capacity(1 << 20, f);
    let mut refs = Reference::default();
    let mut line = Vec::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break;
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
            refs.names.push(name);
            refs.seqs.push(Vec::new());
        } else if !line.is_empty() {
            match refs.seqs.last_mut() {
                Some(s) => s.extend(line.iter().filter(|b| !b.is_ascii_whitespace())),
                None => bail!("sequence data before first FASTA header"),
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
