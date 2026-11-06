#!/usr/bin/env bash
# Run alleleCounter on SNP6 loci, one chromosome per parallel job.
# Writes <out_prefix><chrom>.txt as expected by ascat.getBAFsAndLogRs().
set -euo pipefail

usage() {
    echo "Usage: $0 --bam BAM --ref-dir DIR --out-prefix PREFIX [--chr-prefix true|false] [--min-base-qual 20] [--min-map-qual 35] [--cpus 1] [--chroms '1 ... X']" >&2
    exit 1
}

BAM="" REF_DIR="" OUT_PREFIX="" CHR_PREFIX="false" MIN_BQ=20 MIN_MQ=35 CPUS=1 CHROMS="$(seq -s ' ' 1 22) X"
while [ $# -gt 0 ]; do
    case "$1" in
        --bam)           BAM=$2; shift 2 ;;
        --ref-dir)       REF_DIR=$2; shift 2 ;;
        --out-prefix)    OUT_PREFIX=$2; shift 2 ;;
        --chr-prefix)    CHR_PREFIX=$2; shift 2 ;;
        --min-base-qual) MIN_BQ=$2; shift 2 ;;
        --min-map-qual)  MIN_MQ=$2; shift 2 ;;
        --cpus)          CPUS=$2; shift 2 ;;
        --chroms)        CHROMS=$2; shift 2 ;;
        *)               usage ;;
    esac
done
[ -n "$BAM" ] && [ -n "$REF_DIR" ] && [ -n "$OUT_PREFIX" ] || usage

count_chrom() {
    local c=$1
    local loci="${REF_DIR}/SNP6_loci_chr${c}.txt"
    if [ ! -s "$loci" ]; then
        echo "No SNP6 loci for chromosome ${c} - skipped" >&2
        return 0
    fi
    # match the loci contig names to the BAM ('chr1' vs '1')
    local tmp="loci_chr${c}.tmp"
    if [ "$CHR_PREFIX" = "true" ]; then
        awk 'BEGIN { FS = OFS = "\t" } { if ($1 !~ /^chr/) $1 = "chr" $1; print }' "$loci" > "$tmp"
    else
        awk 'BEGIN { FS = OFS = "\t" } { sub(/^chr/, "", $1); print }' "$loci" > "$tmp"
    fi
    alleleCounter -l "$tmp" -b "$BAM" -o "${OUT_PREFIX}${c}.txt" -m "$MIN_BQ" -q "$MIN_MQ"
    rm -f "$tmp"
}
export -f count_chrom
export BAM REF_DIR OUT_PREFIX CHR_PREFIX MIN_BQ MIN_MQ

# shellcheck disable=SC2086
printf '%s\n' $CHROMS | xargs -P "$CPUS" -I {} bash -c 'count_chrom "$1"' _ {}

n_out=$(find . -maxdepth 1 -name "$(basename "$OUT_PREFIX")*.txt" | wc -l)
echo "alleleCounter finished: ${n_out} chromosome files"
[ "$n_out" -gt 0 ] || { echo "No allele count files produced" >&2; exit 1; }
