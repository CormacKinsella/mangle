#!/usr/bin/env bash
# All pairwise dotplots between the per-product fastas in a mangle output dir.
# Usage: dotplot_all.sh <indir> <outdir> <window>
set -euo pipefail
indir=${1:?indir}; outdir=${2:?outdir}; window=${3:-7}
mkdir -p "$outdir"
# all.fa is the concatenation of every product, so skip it
files=()
for f in "$indir"/*.fa "$indir"/*.fasta "$indir"/*.fna; do
    [[ -e $f && $(basename "$f") != all.* ]] && files+=("$f")
done
n=${#files[@]}
(( n >= 2 )) || { echo "need at least 2 fasta files in $indir" >&2; exit 1; }
for ((i = 0; i < n; i++)); do
    for ((j = i + 1; j < n; j++)); do
        a=$(basename "${files[i]}"); a=${a%.*}
        b=$(basename "${files[j]}"); b=${b%.*}
        dnadotplot --revcompl --svg --window "$window" -o "$outdir/${a}_vs_${b}.svg" \
            -1 "${files[i]}" -2 "${files[j]}"
    done
done
