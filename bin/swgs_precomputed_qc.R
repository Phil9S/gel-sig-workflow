#!/usr/bin/env Rscript
# Per-sample version of swgs-absolutecn v2.0.0 processPrecomputed.R and splitInput.R:
# writes the fit QC row and a one-row sample sheet using Battenberg purity and ploidy.
# Runs in the swgs-absolutecn v2.0.0 container.
#
# rswgsabsolutecn requires 0 < purity <= 1, but Battenberg allows purity up to max_rho
# (1.02 by default), so purity above 1 is capped at 1.

options(scipen = 999)
`%>%` <- dplyr::`%>%`

opts <- optparse::parse_args(
  optparse::OptionParser() |>
    optparse::add_option("--sample", type = "character") |>
    optparse::add_option("--patient", type = "character") |>
    optparse::add_option("--bam", type = "character") |>
    optparse::add_option("--filetype", type = "character", default = "BAM") |>
    optparse::add_option("--bin", type = "double") |>
    optparse::add_option("--genome", type = "character") |>
    optparse::add_option("--ploidy", type = "double") |>
    optparse::add_option("--purity", type = "double") |>
    optparse::add_option("--tp53freq", type = "double", default = NA) |>
    optparse::add_option("--smooth", type = "logical", default = FALSE) |>
    optparse::add_option("--qc_out", type = "character") |>
    optparse::add_option("--sheet_out", type = "character") |>
    optparse::add_option("--bamcheck_out", type = "character")
)

purity <- opts$purity
if (purity > 1) {
  message("[swgs_precomputed_qc] ", opts$sample, ": purity ", purity, " > 1 capped to 1 (rswgsabsolutecn requires purity <= 1)")
  purity <- 1
}
if (purity <= 0 || opts$ploidy <= 1) {
  stop(opts$sample, ": purity (", purity, ") must be > 0 and ploidy (", opts$ploidy, ") > 1 for swgs-absolutecn")
}

rswgsabsolutecn::bamCheck(x = opts$bam, filetype = opts$filetype, outname = opts$bamcheck_out)
if (!file.exists(opts$bamcheck_out)) stop("BAM check failed for ", opts$sample, " (", opts$bam, ")")

samplesheet <- data.frame(PATIENT_ID = opts$patient, SAMPLE_ID = opts$sample, TP53freq = opts$tp53freq,
                          smooth = opts$smooth, file = opts$bam, precPloidy = opts$ploidy,
                          precPurity = purity, stringsAsFactors = FALSE)
write.table(samplesheet, opts$sheet_out, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Fit QC row, as in processPrecomputed.R v2.0.0
bins <- QDNAseqmod::getBinAnnotations(binSize = opts$bin, genome = opts$genome)
# 0.932 ratio between actual usable bins and total usable bins in bin annotation data
nbins_ref_genome <- round(sum(bins@data$use) * 0.932)

fittingColumnNames <- get(utils::data("fittingColumnNames", package = "rswgsabsolutecn", envir = environment()))
preFileCols <- c(fittingColumnNames, "PATIENT_ID", "precPloidy", "precPurity", "TP53freq", "smooth",
                 "rank_clonality", "pl_diff", "new_state_n", "new_state")
exclude <- c("SAMPLE_ID", "PATIENT_ID", "ploidy", "purity", "smooth", "precPloidy", "precPurity", "TP53freq")

samplesheetNew <- samplesheet %>%
  dplyr::mutate(ploidy = precPloidy, purity = precPurity) %>%
  dplyr::relocate(ploidy, purity, .before = "precPloidy") %>%
  dplyr::mutate(!!!setNames(rep(NA, length(preFileCols)), preFileCols[which(!preFileCols %in% exclude)])) %>%
  dplyr::mutate(downsample_depth = rswgsabsolutecn::getDownsampleDepth(ploidy = ploidy, purity = purity,
                                                                         nbins_ref_genome = nbins_ref_genome)) %>%
  dplyr::mutate(use = TRUE) %>%
  dplyr::mutate(notes = "using preprocessed ploidy purity values") %>%
  dplyr::select(dplyr::all_of(preFileCols), use, notes)

write.table(samplesheetNew, opts$qc_out, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
print(samplesheetNew[, c("SAMPLE_ID", "ploidy", "purity", "downsample_depth")])
