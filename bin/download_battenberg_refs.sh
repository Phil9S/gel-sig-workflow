#!/usr/bin/env bash
# Download and unpack the Battenberg v3.0.0 reference files for one genome build into the
# layout conf/genomes.config expects, so battenberg_ref_dir can point straight at --out.
#
#   hg19  ~15.5 GB download. GRCh37 bundle from ORA (Oxford) plus Beagle 5 b37 panels,
#         PLINK maps and jar from the Beagle site (UW), which the ORA bundle lacks.
#   hg38  ~17 GB download. GRCh38 v3.0.0 bundle from Figshare (doi:10.48420/30406441).
#         The bundle ships no impute_info.txt, so one is generated here.
#
# Unpacking needs roughly the same space again; archives are deleted afterwards unless
# --keep-archives is given. Re-running resumes: finished steps are skipped.
set -euo pipefail

usage() {
    echo "Usage: $0 --build hg19|hg38 --out DIR [--keep-archives]" >&2
    exit 1
}

BUILD="" OUT="" KEEP=0
while [ $# -gt 0 ]; do
    case "$1" in
        --build)         BUILD=$2; shift 2 ;;
        --out)           OUT=$2; shift 2 ;;
        --keep-archives) KEEP=1; shift ;;
        *) usage ;;
    esac
done
[ -n "$BUILD" ] && [ -n "$OUT" ] || usage
case "$BUILD" in hg19|hg38) ;; *) usage ;; esac
for tool in curl tar unzip sha1sum; do
    command -v "$tool" >/dev/null || { echo "ERROR: $tool not found" >&2; exit 1; }
done

mkdir -p "$OUT"
cd "$OUT"
mkdir -p .archives

# fetch URL DEST: download with retries; resumes partial files where the server allows it
fetch() {
    local url=$1 dest=$2 rc=0
    if [ -s "$dest" ] || [ -e ".archives/.done.$(basename "$dest")" ]; then
        echo "  have $dest"; return
    fi
    echo "  get  $dest"
    curl -fL -# --retry 5 --retry-delay 15 -C - -o "$dest.part" "$url" || rc=$?
    if [ "$rc" -eq 33 ]; then          # server refused a byte range: start again
        rm -f "$dest.part"
        curl -fL -# --retry 5 --retry-delay 15 -o "$dest.part" "$url"
    elif [ "$rc" -ne 0 ]; then
        echo "ERROR: download failed ($rc): $url" >&2; exit "$rc"
    fi
    mv "$dest.part" "$dest"
}

# unpack ARCHIVE [DIR]: extract once (a marker records completion)
unpack() {
    local a=$1 dir=${2:-.} mark=".archives/.done.$(basename "$1")"
    [ -e "$mark" ] && { echo "  done $a"; return; }
    echo "  unpack $a"
    mkdir -p "$dir"
    case "$a" in
        *.zip)    unzip -q -o "$a" -d "$dir" ;;
        *.tar.gz) tar -xzf "$a" -C "$dir" ;;
    esac
    touch "$mark"
    [ "$KEEP" -eq 1 ] || rm -f "$a"
}

link() {  # link TARGET NAME (inside the current directory)
    [ -e "$1" ] && [ ! -e "$2" ] && ln -s "$1" "$2" || true
}

if [ "$BUILD" = hg19 ]; then
    ORA=https://ora.ox.ac.uk/objects/uuid:2c1fec09-a504-49ab-9ce9-3f17bac531bc/files
    BEAGLE=https://bochet.gcc.biostat.washington.edu/beagle

    echo "[1/3] Battenberg GRCh37 bundle (ORA)"
    while read -r id name; do
        fetch "$ORA/$id" ".archives/$name"
    done <<'EOF'
m684ce995b206d5714ec941e8ed5b30dc battenberg_1000genomesloci2012_v3.tar.gz
r0k225b204 battenberg_impute_1000G_v3.tar.gz
r6969z087n battenberg_wgs_gc_correction_1000g_v3.tar.gz
mfa7e1a13116d8c6b2f5a9b1d0ae79c8b battenberg_wgs_replic_correction_1000g_v3.tar.gz
EOF
    fetch "$ORA/mea5203e1c7d0c3c0bb04d4fda49de56d" probloci_270415.txt.gz
    fetch "$ORA/ma22f8794762c5b8cc473f5f788b8de47" .archives/impute_info.template.txt
    fetch https://raw.githubusercontent.com/Wedge-lab/battenberg/v3.0.0/gcCorrect_chromosome_coordinates_hg19.txt \
          gcCorrect_chromosome_coordinates_hg19.txt
    for a in battenberg_1000genomesloci2012_v3 battenberg_impute_1000G_v3 \
             battenberg_wgs_gc_correction_1000g_v3 battenberg_wgs_replic_correction_1000g_v3; do
        unpack ".archives/$a.tar.gz"
    done

    # the published impute_info.txt holds a placeholder instead of a directory
    sed 's#<path_to_impute_reference_files>#battenberg_impute_v3#g' \
        .archives/impute_info.template.txt > impute_info.txt

    echo "[2/3] Beagle 5 jar, b37 reference panels and PLINK maps"
    mkdir -p beagle/b37
    fetch https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar beagle/beagle.22Jul22.46e.jar
    for c in $(seq 1 22) X; do
        fetch "$BEAGLE/1000_Genomes_phase3_v5a/b37.bref3/chr$c.1kg.phase3.v5a.b37.bref3" \
              "beagle/b37/chr$c.1kg.phase3.v5a.b37.bref3"
    done
    fetch "$BEAGLE/genetic_maps/plink.GRCh37.map.zip" .archives/plink.GRCh37.map.zip
    unpack .archives/plink.GRCh37.map.zip beagle/b37
    # flatten if the zip unpacked into a subfolder
    find beagle/b37 -mindepth 2 -name 'plink.*.map' -exec mv -n {} beagle/b37/ \;

    # Battenberg may ask for chromosome 23 instead of X
    ( cd beagle/b37
      link chrX.1kg.phase3.v5a.b37.bref3 chr23.1kg.phase3.v5a.b37.bref3
      link plink.chrX.GRCh37.map         plink.chr23.GRCh37.map
      link plink.chrX_par1.GRCh37.map    plink.chr23_par1.GRCh37.map
      link plink.chrX_par2.GRCh37.map    plink.chr23_par2.GRCh37.map )
else
    FS=https://ndownloader.figshare.com/files

    echo "[1/3] Battenberg GRCh38 v3.0.0 bundle (Figshare)"
    while read -r id name; do
        fetch "$FS/$id" "$name"
    done <<'EOF'
58912789 sha1sums.txt
58912798 README.txt
58912792 chromosome_coordinates_hg38.txt
58912795 chromosome_coordinates_hg38_chr.txt
58912801 probloci.zip
58912804 imputation.zip
58912807 RT_correction_hg38.zip
58912810 1000G_loci_hg38.zip
58912813 GC_correction_hg38.zip
58912840 shapeit2.zip
58912852 beagle.zip
EOF
    echo "  verifying checksums"
    # skip archives already unpacked (and deleted) by an earlier run
    to_check=$(while read -r sum name; do
        [ -e ".archives/.done.$name" ] || echo "$sum  $name"
    done < sha1sums.txt)
    [ -z "$to_check" ] || echo "$to_check" | sha1sum -c --quiet -
    for a in probloci imputation RT_correction_hg38 1000G_loci_hg38 GC_correction_hg38 shapeit2 beagle; do
        unpack "$a.zip"
    done

    echo "[2/3] Generating imputation/impute_info.txt"
    # columns: chrom, legend, genetic map, haplotypes, start, end, is_par (1 = diploid in males)
    S=shapeit2/ALL.v1a.shapeit2_integrated_chr
    M=imputation/genetic_map_chr
    {
        awk 'NR > 1 && $1 <= 22 {print $1, $5}' chromosome_coordinates_hg38.txt |
        while read -r c end; do
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$c" "$S$c.GRCh38.20181129.phased.legend" \
                "${M}${c}_combined_b38.txt" "$S$c.GRCh38.20181129.phased.hap.gz" 0 "$end" 1
        done
        # GRCh38 pseudo-autosomal boundaries
        for r in "PAR1 0 2781479 1" "nonPAR 2781480 155701382 0" "PAR2 155701383 156040895 1"; do
            set -- $r
            printf 'X\t%s\t%s\t%s\t%s\t%s\t%s\n' "${S}X_$1.GRCh38.20181129.phased.legend" \
                "${M}X_$1_combined_b38.txt" "${S}X_$1.GRCh38.20181129.phased.hap.gz" "$2" "$3" "$4"
        done
    } > imputation/impute_info.txt
fi

echo "[3/3] Checking files referenced by impute_info.txt"
INFO=$([ "$BUILD" = hg19 ] && echo impute_info.txt || echo imputation/impute_info.txt)
missing=$(cut -f2-4 "$INFO" | tr '\t' '\n' | while read -r f; do [ -e "$f" ] || echo "$f"; done)
if [ -n "$missing" ]; then
    echo "ERROR: files listed in $INFO are missing:" >&2
    echo "$missing" | head >&2
    exit 1
fi

echo
echo "Done. In your config set:"
echo "    params.genomes.$BUILD.battenberg_ref_dir = '$(pwd -P)'"
