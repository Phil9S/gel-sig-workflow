# ASCAT v3.2.0 and alleleCounter, built from the locked pixi environment, with ASCAT's WGS
# 1000 Genomes references (loci, alleles, GC, RT) for hg19 and hg38 in /ref/<build>/, plus the
# SNP6 positions: /ref/hg19/SNPpos.txt and /ref/hg38/SNPpos_hg38.txt (lifted over at build).
# Tag <ASCAT version>-<image revision>; bump the revision when only the image changes.
# docker build -t phil9s/ascat:3.2.0-1 containers/ascat
# docker push phil9s/ascat:3.2.0-1
# Run `pixi lock` after editing pixi.toml.

FROM ghcr.io/prefix-dev/pixi:0.82.1 AS build
ARG ASCAT_REF=v3.2.0
WORKDIR /app
COPY pixi.toml pixi.lock ./
# post-link scripts install Bioconductor data packages (GenomeInfoDbData)
RUN pixi install --locked --run-post-link-scripts
# dependencies = FALSE: R deps come from pixi.lock, not CRAN. install_github only warns when
# the install fails, so check the package loads.
RUN export TAR="tar --no-same-permissions --no-same-owner" TMPDIR=/tmp && \
    pixi run Rscript -e "remotes::install_github(repo = 'VanLoo-lab/ascat/ASCAT@${ASCAT_REF}', quiet = TRUE, upgrade = FALSE, dependencies = FALSE); if (!requireNamespace('ASCAT', quietly = TRUE)) quit(status = 1)"
# r-base pulls in GCC/gfortran/binutils, sanitizer runtimes and a glibc sysroot to compile
# packages; Python only runs the data-package post-link script. ASCAT is pure R and already
# installed, so drop them, along with static libraries, headers and documentation.
RUN cd /app/.pixi/envs/default && \
    rm -rf libexec/gcc lib/gcc x86_64-conda-linux-gnu bin/x86_64-conda-linux-gnu-* \
           lib/libasan.* lib/libtsan.* lib/libhwasan.* lib/liblsan.* lib/libubsan.* lib/*.a \
           lib/python3* bin/python* \
           include share/doc share/man share/info share/gtk-doc lib/R/doc/manual \
           lib/R/library/*/doc lib/R/library/*/html

# ASCAT WGS reference files and SNP6 positions: https://zenodo.org/records/14008443
FROM debian:12-slim AS refs
ARG LIFTOVER_VERSION=482
ARG CHAIN_MD5=35887f73fe5e2231656504d1f6430900
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates curl unzip && \
    rm -rf /var/lib/apt/lists/*
COPY --from=ghcr.io/prefix-dev/pixi:0.82.1 /usr/local/bin/pixi /usr/local/bin/pixi
WORKDIR /dl
COPY ref_files.md5 lift_snppos.sh ./
# <type>_WGS_<build>.zip unpacks to /ref/<build>/; SNPpos.txt (SNP6 positions) is hg19
RUN while read -r sum f; do \
        curl -fsSL --retry 5 --retry-delay 10 -o "$f" "https://zenodo.org/records/14008443/files/$f?download=1" || exit 1; \
    done < ref_files.md5 && \
    md5sum -c ref_files.md5 && \
    for f in *.zip; do \
        case "$f" in SNPpos.txt.zip) b=hg19 ;; *) b=${f%.zip}; b=${b##*_} ;; esac; \
        mkdir -p "/ref/$b" && unzip -q -o "$f" -d "/ref/$b" || exit 1; \
    done && \
    rm -f *.zip
# hg38 SNP6 positions: no GRCh38 file is published, so lift SNPpos.txt over with the UCSC chain
RUN curl -fsSL --retry 5 --retry-delay 10 -o hg19ToHg38.over.chain.gz \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz && \
    echo "${CHAIN_MD5}  hg19ToHg38.over.chain.gz" | md5sum -c && \
    pixi exec -c conda-forge -c bioconda -s "ucsc-liftover==${LIFTOVER_VERSION}" \
        bash lift_snppos.sh /ref/hg19/SNPpos.txt hg19ToHg38.over.chain.gz /ref/hg38/SNPpos_hg38.txt && \
    chmod -R a+rX /ref

FROM debian:12-slim AS production
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8
# reference data first: it changes least often, so its layer is reused across rebuilds
COPY --from=refs /ref /ref
WORKDIR /app
COPY pixi.toml pixi.lock ./
# The env prefix must match the build stage
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
ENV PATH=/app/.pixi/envs/default/bin:$PATH \
    CONDA_PREFIX=/app/.pixi/envs/default

# smoke test: ASCAT loads, cairo PNG output works, alleleCounter runs, references are present
RUN Rscript -e "suppressPackageStartupMessages(library(ASCAT)); cat('ASCAT', as.character(packageVersion('ASCAT')), '\n'); png(tempfile(fileext = '.png'), type = 'cairo'); plot(1); invisible(dev.off())" && \
    alleleCounter --version && \
    for b in hg19 hg38; do \
        snppos=SNPpos.txt; [ "$b" = hg19 ] || snppos=SNPpos_$b.txt; \
        for f in "G1000_loci_${b}_chr1.txt" "G1000_alleles_${b}_chrX.txt" "GC_G1000_${b}.txt" "RT_G1000_${b}.txt" "$snppos"; do \
            test -s "/ref/$b/$f" || { echo "missing /ref/$b/$f" >&2; exit 1; }; \
        done; \
    done

LABEL org.opencontainers.image.title="ascat" \
      org.opencontainers.image.version="3.2.0-1" \
      org.opencontainers.image.source="https://github.com/VanLoo-lab/ascat" \
      org.opencontainers.image.description="ASCAT v3.2.0 + alleleCounter with WGS G1000 references and SNP6 positions for hg19/hg38 in /ref/<build>/"
