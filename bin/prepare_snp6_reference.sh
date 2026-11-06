#!/usr/bin/env bash
# Build ASCAT HTS loci, alleles, GC and RT reference files restricted to SNP6.0 positions.
# Positions are intersected with ASCAT's 1000 Genomes references so all files stay consistent.
# The positions file is matched on chr/pos header names, else the first two columns
# (columns 2-3 when column 1 is a probe ID).
set -euo pipefail

usage() {
    echo "Usage: $0 --positions FILE --loci-prefix PREFIX --alleles-prefix PREFIX --gc FILE [--rt FILE] --out DIR [--chroms '1 2 ... X']" >&2
    exit 1
}

POSITIONS="" LOCI="" ALLELES="" GC="" RT="" OUT="" CHROMS="$(seq -s ' ' 1 22) X"
while [ $# -gt 0 ]; do
    case "$1" in
        --positions)      POSITIONS=$2; shift 2 ;;
        --loci-prefix)    LOCI=$2; shift 2 ;;
        --alleles-prefix) ALLELES=$2; shift 2 ;;
        --gc)             GC=$2; shift 2 ;;
        --rt)             RT=$2; shift 2 ;;
        --out)            OUT=$2; shift 2 ;;
        --chroms)         CHROMS=$2; shift 2 ;;
        *)                usage ;;
    esac
done
[ -n "$POSITIONS" ] && [ -n "$LOCI" ] && [ -n "$ALLELES" ] && [ -n "$GC" ] && [ -n "$OUT" ] || usage
[ -s "$GC" ] || { echo "GC correction file not found: $GC" >&2; exit 1; }
mkdir -p "$OUT"
KEYS="$OUT/snp6_keys.txt"

# SNP6 positions as <chr>_<pos> keys
zcat -f "$POSITIONS" | awk -F'\t' '
    function norm(c) { gsub(/"/, "", c); sub(/^chr/, "", c); if (c == "23") c = "X"; if (c == "24") c = "Y"; return c }
    function is_chrom(c) { c = norm(c); return c ~ /^([0-9]+|X|Y)$/ }
    { sub(/\r$/, "") }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            h = tolower($i); gsub(/"/, "", h)
            if (!ci && (h == "chr" || h == "chrom" || h == "chromosome")) ci = i
            if (!pi && (h == "pos" || h == "position" || h == "physical_position" || h == "start")) pi = i
        }
        if (ci && pi) next
        if (is_chrom($1) && $2 ~ /^[0-9]+$/) { ci = 1; pi = 2 }
        else if (is_chrom($2) && $3 ~ /^[0-9]+$/) { ci = 2; pi = 3 }
        else { print "Cannot identify chromosome/position columns in SNP6 positions file" > "/dev/stderr"; exit 2 }
    }
    { p = $pi; gsub(/"/, "", p); if (p ~ /^[0-9]+$/ && is_chrom($ci)) print norm($ci) "_" p }
' | sort -u > "$KEYS"
N_POSITIONS=$(wc -l < "$KEYS")
echo "SNP6 positions: $N_POSITIONS"
[ "$N_POSITIONS" -gt 0 ] || { echo "No positions read from $POSITIONS" >&2; exit 1; }

# Loci and alleles per chromosome
printf 'chromosome\tsnp6_positions\tloci_retained\n' > "$OUT/snp6_reference_summary.tsv"
TOTAL=0
for c in $CHROMS; do
    loci_in="${LOCI}${c}.txt"
    alleles_in="${ALLELES}${c}.txt"
    loci_out="$OUT/SNP6_loci_chr${c}.txt"
    alleles_out="$OUT/SNP6_alleles_chr${c}.txt"
    n_pos=$(awk -v c="$c" 'index($0, c "_") == 1' "$KEYS" | wc -l)
    if [ ! -s "$loci_in" ] || [ ! -s "$alleles_in" ]; then
        echo "WARNING: missing $loci_in or $alleles_in - chromosome $c skipped" >&2
        printf '%s\t%s\t0\n' "$c" "$n_pos" >> "$OUT/snp6_reference_summary.tsv"
        continue
    fi
    awk -v c="$c" 'BEGIN { FS = OFS = "\t" }
        FILENAME == ARGV[1] { if (index($0, c "_") == 1) k[$0] = 1; next }
        { chr = $1; sub(/^chr/, "", chr); if ((chr "_" $2) in k) print }' "$KEYS" "$loci_in" > "$loci_out"
    awk 'BEGIN { FS = OFS = "\t" }
        FILENAME == ARGV[1] { k[$2] = 1; next }
        FNR == 1 { print; next }
        ($1 in k)' "$loci_out" "$alleles_in" > "$alleles_out"
    n_loci=$(wc -l < "$loci_out")
    TOTAL=$((TOTAL + n_loci))
    printf '%s\t%s\t%s\n' "$c" "$n_pos" "$n_loci" >> "$OUT/snp6_reference_summary.tsv"
done
echo "SNP6 loci retained in the 1000 Genomes allele reference: $TOTAL"
[ "$TOTAL" -ge 1000 ] || { echo "Too few SNP6 loci retained ($TOTAL) - check chromosome naming and build of the inputs" >&2; exit 1; }

# GC and RT correction files
subset_correction() {
    local in=$1 out=$2
    cat "$OUT"/SNP6_loci_chr*.txt | awk 'BEGIN { OFS = "\t" } { sub(/^chr/, "", $1); print $1 "_" $2 }' > "$OUT/retained_keys.txt"
    zcat -f "$in" | awk 'BEGIN { FS = OFS = "\t" }
        function norm(c) { gsub(/"/, "", c); sub(/^chr/, "", c); return c }
        FILENAME == ARGV[1] { k[$0] = 1; next }
        { sub(/\r$/, "") }
        FNR == 1 { print; next }
        { key = norm($2) "_" $3; if (key in k) { $1 = key; print } }' "$OUT/retained_keys.txt" - > "$out"
    rm -f "$OUT/retained_keys.txt"
    echo "$(basename "$out"): $(( $(wc -l < "$out") - 1 )) loci"
}
subset_correction "$GC" "$OUT/GC_SNP6.txt"
if [ -n "$RT" ] && [ -s "$RT" ]; then
    subset_correction "$RT" "$OUT/RT_SNP6.txt"
else
    echo "No replication timing file - RT correction will be skipped"
fi
rm -f "$KEYS"
