#!/usr/bin/env Rscript
# Extract Battenberg purity, ploidy, psi and goodness of fit into a one-row table.
# ploidy is the segment-weighted ploidy of the final profile (used for sWGS absolute CN);
# psi is the fit parameter ASCAT expects as psi_manual.

library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample", type = "character", default = NA),
  make_option("--tumour_name", type = "character", default = NA),
  make_option("--normal_name", type = "character", default = NA),
  make_option("--cancer", type = "character", default = NA),
  make_option("--build", type = "character", default = NA),
  make_option("--sex", type = "character", default = NA),
  make_option("--purity_ploidy", type = "character", default = NA, help = "Battenberg _purity_ploidy.txt"),
  make_option("--rho_psi", type = "character", default = NA, help = "Battenberg _rho_and_psi.txt"),
  make_option("--out", type = "character", help = "output TSV")
)))
if (is.null(opt$out)) stop("--out is required")
get <- function(key) opt[[key]]

num_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  suppressWarnings(as.numeric(x[1]))
}

pp <- NULL
pp_file <- get("purity_ploidy")
if (!is.na(pp_file) && file.exists(pp_file) && file.info(pp_file)$size > 0) {
  pp <- read.delim(pp_file, stringsAsFactors = FALSE)
  colnames(pp)[colnames(pp) == "cellularity"] <- "purity"
}

rp <- NULL
rp_file <- get("rho_psi")
if (!is.na(rp_file) && file.exists(rp_file) && file.info(rp_file)$size > 0) {
  # header has one field fewer than the rows: the first column holds the row names
  rp <- read.table(rp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}
from_rp <- function(col, row = "FRAC_GENOME") {
  if (is.null(rp) || !(row %in% rownames(rp)) || !(col %in% colnames(rp))) return(NA_real_)
  num_or_na(rp[row, col])
}

purity <- if (!is.null(pp) && "purity" %in% colnames(pp)) num_or_na(pp$purity) else NA_real_
ploidy <- if (!is.null(pp) && "ploidy" %in% colnames(pp)) num_or_na(pp$ploidy) else NA_real_
psi    <- if (!is.null(pp) && "psi" %in% colnames(pp)) num_or_na(pp$psi) else NA_real_

if (is.na(purity)) purity <- from_rp("rho")
if (is.na(ploidy)) ploidy <- from_rp("ploidy")
if (is.na(psi))    psi    <- from_rp("psi")
goodness <- from_rp("distance")

out <- data.frame(
  sample          = get("sample"),
  tumour_name     = get("tumour_name"),
  normal_name     = get("normal_name"),
  cancer          = get("cancer"),
  build           = get("build"),
  sex             = get("sex"),
  purity          = round(purity, 4),
  ploidy          = round(ploidy, 4),
  psi             = round(psi, 4),
  goodness_of_fit = round(goodness, 4),
  solution        = "FRAC_GENOME",
  stringsAsFactors = FALSE
)
if (is.na(out$purity) || is.na(out$ploidy)) {
  message("WARNING: no purity/ploidy could be extracted for sample ", out$sample)
}
write.table(out, get("out"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
print(out)
