#!/usr/bin/env bash
# Lift ASCAT's SNP6 positions (SNPpos.txt: probe ID, Chr, Position; hg19) to hg38 with UCSC
# liftOver. Keeps loci that map to a single base on the same chromosome (1-22, X) and writes
# the same three-column layout without a 'chr' prefix.
# Usage: lift_snppos.sh SNPpos.txt hg19ToHg38.over.chain.gz OUT.txt
set -euo pipefail
[ $# -eq 3 ] || { echo "Usage: $0 SNPpos.txt CHAIN OUT.txt" >&2; exit 1; }
in=$1 chain=$2 out=$3
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# 0-based BED, 'chr' names for the UCSC chain; the name keeps probe ID and source chromosome
awk -F'\t' 'BEGIN { OFS = "\t" }
    { sub(/\r$/, "") }
    NR > 1 && $3 ~ /^[0-9]+$/ { print "chr" $2, $3 - 1, $3, $1 "|" $2 }' "$in" > "$tmp/in.bed"
liftOver "$tmp/in.bed" "$chain" "$tmp/out.bed" "$tmp/unmapped.bed"

printf '\tChr\tPosition\n' > "$out"
awk -F'\t' 'BEGIN { OFS = "\t" }
    { c = $1; sub(/^chr/, "", c); split($4, a, "|")
      if (c == a[2] && c ~ /^([0-9]+|X)$/ && $3 - $2 == 1) print a[1], c, $3 }' "$tmp/out.bed" |
    sort -t$'\t' -k2,2V -k3,3n >> "$out"

n_in=$(wc -l < "$tmp/in.bed")
n_out=$(( $(wc -l < "$out") - 1 ))
echo "SNP6 positions lifted hg19 -> hg38: $n_out of $n_in kept ($(( n_in - n_out )) unmapped or moved chromosome)"
# a wrong chain or input build drops far more than the usual ~1%
[ $(( n_out * 100 )) -ge $(( n_in * 95 )) ] || { echo "Fewer than 95% of positions lifted over" >&2; exit 1; }
