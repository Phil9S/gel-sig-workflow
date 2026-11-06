#!/usr/bin/env Rscript
# Join the cohort tables (Battenberg, VAF evaluation, ASCAT SNP6, swgs-absolutecn fits)
# into one row per sample. Only the Battenberg table is required.

library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--tables_dir", type = "character", help = "directory holding the cohort tables"),
  make_option("--out", type = "character", help = "output TSV")
)))
if (is.null(opt$tables_dir) || is.null(opt$out)) stop("--tables_dir and --out are required")
options(scipen = 999)

read_table <- function(name) {
  f <- file.path(opt$tables_dir, name)
  if (!file.exists(f) || file.info(f)$size == 0) return(NULL)
  # read as text so IDs such as "T" or "NA" are not converted, then type the other columns
  df <- read.delim(f, check.names = FALSE, colClasses = "character", na.strings = character(0))
  ids <- intersect(c("sample", "tumour_name", "normal_name", "cancer", "SAMPLE_ID", "PATIENT_ID"), colnames(df))
  for (col in setdiff(colnames(df), ids)) df[[col]] <- type.convert(df[[col]], as.is = TRUE, na.strings = "NA")
  df
}
pick <- function(df, cols, prefix = "") {
  cols <- intersect(cols, colnames(df))
  out <- df[, c("sample", cols), drop = FALSE]
  colnames(out)[-1] <- paste0(prefix, cols)
  out
}
merge_in <- function(base, add) if (is.null(add)) base else merge(base, add, by = "sample", all.x = TRUE, sort = FALSE)

bb <- read_table("battenberg_purity_ploidy.tsv")
if (is.null(bb)) stop("battenberg_purity_ploidy.tsv not found")
summary <- bb
colnames(summary)[colnames(summary) %in% c("purity", "ploidy", "psi", "goodness_of_fit")] <-
  paste0("battenberg_", colnames(summary)[colnames(summary) %in% c("purity", "ploidy", "psi", "goodness_of_fit")])
summary$solution <- NULL

vaf <- read_table("vaf_evaluation_metrics.tsv")
if (!is.null(vaf)) {
  summary <- merge_in(summary, pick(vaf, c("qc_flag", "qc_reason", "n_snvs_clonal_segments", "ccf_clonal_peak",
                                           "frac_ccf_near_1", "frac_exceed_max", "snv_implied_purity",
                                           "purity_delta", "frac_multiplicity_gt1"), "vaf_"))
}

ascat <- read_table("ascat_snp6_metrics.tsv")
if (!is.null(ascat)) {
  summary <- merge_in(summary, pick(ascat, c("ascat_fixed_ploidy", "ascat_fixed_goodness_of_fit", "ascat_fixed_failed",
                                             "ascat_unconstrained_purity", "ascat_unconstrained_ploidy",
                                             "ascat_unconstrained_goodness_of_fit", "ascat_qc_n_het_SNP",
                                             "ascat_qc_WGD", "ascat_qc_GI", "ascat_qc_LOH")))
}

fits_file <- list.files(opt$tables_dir, pattern = "_ds_abs_fits\\.tsv$")
fits <- if (length(fits_file) == 1) read_table(fits_file) else NULL
if (!is.null(fits)) {
  # .post columns are the fit of the downsampled profile; purity is capped at 1
  fits$sample <- as.character(fits$SAMPLE_ID)
  summary <- merge_in(summary, pick(fits, c("purity", "ploidy", "downsample_depth", "total.reads", "used.reads",
                                            "segments.post", "clonality.post", "rmse.post", "segvariance.post",
                                            "TP53cn.post", "expected_TP53_AF.post"), "swgs_"))
}

summary <- summary[order(summary$cancer, summary$sample), , drop = FALSE]
write.table(summary, opt$out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Summary written for", nrow(summary), "samples\n")
