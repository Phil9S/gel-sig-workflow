#!/usr/bin/env Rscript
# Run ASCAT at SNP6.0 loci with purity (rho) and ploidy (psi) fixed to the Battenberg fit,
# optionally alongside ASCAT's unconstrained fit for comparison.
# Writes segments, a metrics table with ascat.metrics() QC, an RData file and plots.

suppressPackageStartupMessages(library(ASCAT))
options(bitmapType = "cairo", stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (i == length(args) || startsWith(args[i + 1], "--")) {
      out[[key]] <- "TRUE"; i <- i + 1
    } else {
      out[[key]] <- args[i + 1]; i <- i + 2
    }
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
flag <- function(x) !is.null(x) && toupper(x) %in% c("TRUE", "T", "1", "YES")

SAMPLE      <- opt$sample
TUMOUR_NAME <- opt$tumour_name
NORMAL_NAME <- opt$normal_name
SEX         <- opt$sex
GENOME      <- opt$build
PURITY      <- as.numeric(opt$purity)
PSI         <- as.numeric(opt$psi)
BB_PLOIDY   <- as.numeric(opt$ploidy)
GAMMA       <- as.numeric(opt$gamma)
PENALTY     <- as.numeric(opt$penalty)
MIN_COUNTS  <- as.numeric(opt$min_counts)
SEED        <- as.integer(opt$seed)
GC_FILE     <- opt$gc_file
RT_FILE     <- if (!is.null(opt$rt_file) && file.exists(opt$rt_file)) opt$rt_file else NULL

cat("[run_ascat_snp6] ASCAT", as.character(packageVersion("ASCAT")), "| sample", SAMPLE,
    "| fixed purity", PURITY, "psi", PSI, "\n")

chrom_names <- c(1:22, "X")
available <- chrom_names[file.exists(paste0(opt$tumour_counts_prefix, chrom_names, ".txt")) &
                         file.exists(paste0(opt$normal_counts_prefix, chrom_names, ".txt")) &
                         file.exists(paste0(opt$alleles_prefix, chrom_names, ".txt"))]
if (length(available) == 0) stop("No chromosomes with tumour, normal and allele files")
if (!identical(available, chrom_names)) cat("[run_ascat_snp6] chromosomes used:", paste(available, collapse = ","), "\n")

tumourLogR_file  <- paste0(TUMOUR_NAME, "_tumour_LogR.txt")
tumourBAF_file   <- paste0(TUMOUR_NAME, "_tumour_BAF.txt")
normalLogR_file  <- paste0(TUMOUR_NAME, "_germline_LogR.txt")
normalBAF_file   <- paste0(TUMOUR_NAME, "_germline_BAF.txt")

ascat.getBAFsAndLogRs(samplename = TUMOUR_NAME,
                      tumourAlleleCountsFile.prefix = opt$tumour_counts_prefix,
                      normalAlleleCountsFile.prefix = opt$normal_counts_prefix,
                      tumourLogR_file = tumourLogR_file,
                      tumourBAF_file = tumourBAF_file,
                      normalLogR_file = normalLogR_file,
                      normalBAF_file = normalBAF_file,
                      alleles.prefix = opt$alleles_prefix,
                      gender = SEX,
                      genomeVersion = GENOME,
                      chrom_names = available,
                      minCounts = MIN_COUNTS,
                      seed = SEED)
ascat.synchroniseFiles(samplename = TUMOUR_NAME,
                       tumourLogR_file = tumourLogR_file,
                       tumourBAF_file = tumourBAF_file,
                       normalLogR_file = normalLogR_file,
                       normalBAF_file = normalBAF_file)

ascat.bc <- ascat.loadData(Tumor_LogR_file = tumourLogR_file,
                           Tumor_BAF_file = tumourBAF_file,
                           Germline_LogR_file = normalLogR_file,
                           Germline_BAF_file = normalBAF_file,
                           gender = SEX,
                           genomeVersion = GENOME)
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc <- ascat.correctLogR(ascat.bc, GCcontentfile = GC_FILE, replictimingfile = RT_FILE)
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc <- ascat.aspcf(ascat.bc, penalty = PENALTY, seed = SEED)
ascat.plotSegmentedData(ascat.bc)

# Fit with purity and psi fixed to Battenberg
ascat.output <- ascat.runAscat(ascat.bc, gamma = GAMMA, rho_manual = PURITY, psi_manual = PSI,
                               img.prefix = "Battenberg_fixed_", write_segments = FALSE)

write_segments <- function(fit, suffix) {
  if (length(fit$failedarrays) > 0 || is.null(fit$segments)) {
    seg <- data.frame(sample = SAMPLE, chr = NA, startpos = NA, endpos = NA, nMajor = NA, nMinor = NA)
    seg_raw <- cbind(seg, nAraw = NA, nBraw = NA)
  } else {
    seg <- fit$segments
    seg$sample <- SAMPLE
    seg_raw <- fit$segments_raw
    seg_raw$sample <- SAMPLE
  }
  write.table(seg, paste0(SAMPLE, suffix, ".segments.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(seg_raw, paste0(SAMPLE, suffix, ".segments_raw.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}
write_segments(ascat.output, ".ascat")

first_or_na <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x[1])
metrics <- data.frame(
  sample              = SAMPLE,
  tumour_name         = TUMOUR_NAME,
  battenberg_purity   = PURITY,
  battenberg_psi      = PSI,
  battenberg_ploidy   = BB_PLOIDY,
  ascat_fixed_purity  = first_or_na(ascat.output$purity),
  ascat_fixed_ploidy  = first_or_na(ascat.output$ploidy),
  ascat_fixed_psi     = first_or_na(ascat.output$psi),
  ascat_fixed_goodness_of_fit = first_or_na(ascat.output$goodnessOfFit),
  ascat_fixed_failed  = length(ascat.output$failedarrays) > 0
)

# ascat.metrics() columns (ASCAT 3.1.x/3.2.x); always written so per-sample tables share a header
qc_cols <- c("sex", "tumour_mapd", "normal_mapd", "GC_correction_before", "GC_correction_after",
             "RT_correction_before", "RT_correction_after", "n_het_SNP", "n_segs_logR", "n_segs_BAF",
             "n_segs_logRBAF_diff", "frac_homo", "size_intermediate_segments", "size_odd_segments", "n_segs",
             "segs_size", "homdel_segs", "homdel_largest", "homdel_size", "homdel_fraction", "LOH",
             "mode_minA", "mode_majA", "WGD", "GI")
qc <- tryCatch(ascat.metrics(ascat.bc, ascat.output), error = function(e) {
  message("ascat.metrics failed: ", conditionMessage(e)); NULL
})
for (col in qc_cols) {
  value <- if (!is.null(qc) && nrow(qc) > 0 && col %in% colnames(qc)) qc[1, col] else NA
  metrics[[paste0("ascat_qc_", col)]] <- value
}

# ASCAT's unconstrained fit, for comparison
ascat.output.unconstrained <- NULL
if (flag(opt$run_unconstrained)) {
  ascat.output.unconstrained <- tryCatch(
    ascat.runAscat(ascat.bc, gamma = GAMMA, img.prefix = "Unconstrained_", write_segments = FALSE),
    error = function(e) { message("Unconstrained ASCAT fit failed: ", conditionMessage(e)); NULL })
  if (!is.null(ascat.output.unconstrained)) write_segments(ascat.output.unconstrained, ".ascat_unconstrained")
}
metrics$ascat_unconstrained_purity <- first_or_na(ascat.output.unconstrained$purity)
metrics$ascat_unconstrained_ploidy <- first_or_na(ascat.output.unconstrained$ploidy)
metrics$ascat_unconstrained_goodness_of_fit <- first_or_na(ascat.output.unconstrained$goodnessOfFit)

write.table(metrics, paste0(SAMPLE, ".ascat.metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
save(ascat.bc, ascat.output, ascat.output.unconstrained, metrics, file = paste0(SAMPLE, ".ascat.RData"))
print(t(metrics))
sessionInfo()
