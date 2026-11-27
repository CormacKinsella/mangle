#!/usr/bin/env bash
# Full variant list per product from a mangle output dir, whatever form products.tsv uses (ids or row ranges).
# Usage: product_variants.sh <outdir> [product]
# With a product: one variant id per line. Without: a two-column table, product<TAB>variant, for every product.
set -euo pipefail
dir=${1:?outdir}; want=${2:-}
awk -F'\t' -v want="$want" '
    FNR == NR { if (FNR > 1) id[FNR - 1] = $1; next }   # resolved_catalogue.tsv: row number -> id
    FNR == 1 { ranges = ($4 == "variant_rows"); next }
    want != "" && $1 != want { next }
    $4 == "." { next }
    {
        n = split($4, part, ",")
        for (i = 1; i <= n; i++) {
            if (!ranges) { emit($1, part[i]); continue }
            split(part[i], be, "-"); lo = be[1]; hi = (be[2] == "" ? be[1] : be[2])
            for (r = lo; r <= hi; r++) emit($1, id[r])
        }
    }
    function emit(p, v) { print (want == "" ? p "\t" v : v) }
' "$dir/resolved_catalogue.tsv" "$dir/products.tsv"
