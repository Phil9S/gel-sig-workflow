# Battenberg v3.0.0 image built from the locked pixi environment; reference data and the
# Beagle jar come from the Battenberg reference bundle, not the image.
# docker build -t phil9s/battenberg:3.0.0 containers/battenberg
# docker push phil9s/battenberg:3.0.0
# Run `pixi lock` after editing pixi.toml.
FROM ghcr.io/prefix-dev/pixi:0.82.1 AS build
ARG BATTENBERG_REF=v3.0.0
ARG ASCAT_REF=v3.2.0
ARG COPYNUMBER_REF=c00fce95dcce7b1d56b902459945e4bc3f55b1b4
WORKDIR /app
COPY pixi.* .
RUN pixi install --locked --run-post-link-scripts
# dependencies = FALSE: R deps come from pixi.lock, not CRAN
RUN export TAR="tar --no-same-permissions --no-same-owner" && \
	export TMPDIR=/tmp && \
	pixi run Rscript -e "remotes::install_github(repo = 'igordot/copynumber@${COPYNUMBER_REF}', quiet = TRUE, upgrade = FALSE, force = FALSE, dependencies = FALSE)"
RUN export TAR="tar --no-same-permissions --no-same-owner" && \
	export TMPDIR=/tmp && \
	pixi run Rscript -e "remotes::install_github(repo = 'VanLoo-lab/ascat/ASCAT@${ASCAT_REF}', quiet = TRUE, upgrade = FALSE, force = FALSE, dependencies = FALSE)"
RUN export TAR="tar --no-same-permissions --no-same-owner" && \
	export TMPDIR=/tmp && \
	pixi run Rscript -e "remotes::install_github(repo = 'Wedge-lab/battenberg@${BATTENBERG_REF}', quiet = TRUE, upgrade = FALSE, force = FALSE, dependencies = FALSE)"
RUN pixi run Rscript -e "for (p in c('copynumber', 'ASCAT', 'Battenberg')) { \
	  suppressPackageStartupMessages(library(p, character.only = TRUE)); \
	  cat(p, as.character(packageVersion(p)), '\n') }"

FROM debian:12-slim AS production
RUN TZ=Etc/UTC && \
ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install -y --no-install-recommends locales && \
    rm -rf /var/lib/apt/lists/*

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV LANGUAGE=C.UTF-8
WORKDIR /app
COPY pixi.toml pixi.lock ./
# The env prefix must match the build stage
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default

ENV PATH=/app/.pixi/envs/default/bin:$PATH
ENV CONDA_PREFIX=/app/.pixi/envs/default

LABEL org.opencontainers.image.title="battenberg" \
      org.opencontainers.image.version="3.0.0" \
      org.opencontainers.image.source="https://github.com/Wedge-lab/battenberg"
