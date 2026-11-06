#!/usr/bin/env Rscript
# Evaluate a Battenberg copy number fit with somatic SNV VAFs on clonal segments.
# Writes per-sample metrics, per-state peaks, per-SNV values and QC plots.
#
# A clonal SNV on m of N = nMaj + nMin copies at purity p has expected
#   VAF(m) = m * p / (p * N + 2 * (1 - p)),  m = 1..nMaj
# so with correct purity and copy number, CCF clusters at 1. A low fraction of
# multiplicity > 1 SNVs on nMaj >= 2 segments can indicate a spurious genome doubling.
#
# VCF counts: Strelka2 tier-1 AU/CU/GU/TU, else AD, else AF * DP.

suppressPackageStartupMessages({
  library(optparse)
  library(vcfR)
  library(dplyr)
  library(ggplot2)
})
options(stringsAsFactors = FALSE, warn = 1)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample", type = "character", help = "sample ID"),
  make_option("--tumour_name", type = "character", default = "", help = "tumour sample name in the VCF"),
  make_option("--normal_name", type = "character", default = "", help = "normal sample name in the VCF"),
  make_option("--copynumber", type = "character", help = "Battenberg _copynumber.txt"),
  make_option("--purity_ploidy", type = "character", help = "purity/ploidy table from extract_purity_ploidy.R"),
  make_option("--vcf", type = "character", help = "somatic SNV VCF (plain or gzipped)"),
  make_option("--min_depth", type = "double", default = 20),
  make_option("--min_alt_reads", type = "double", default = 3),
  make_option("--pass_only", type = "logical", default = TRUE),
  make_option("--autosomes_only", type = "logical", default = TRUE),
  make_option("--max_total_cn", type = "double", default = 8),
  make_option("--min_snvs", type = "double", default = 100),
  make_option("--min_snvs_per_state", type = "double", default = 50),
  make_option("--ccf_tolerance", type = "double", default = 0.1),
  make_option("--max_exceed_fraction", type = "double", default = 0.05),
  make_option("--out_prefix", type = "character", default = NULL, help = "default <sample>.vaf_eval")
)))
missing_args <- setdiff(c("sample", "copynumber", "purity_ploidy", "vcf"), names(opt))
if (length(missing_args)) stop("Missing required arguments: ", paste0("--", missing_args, collapse = ", "))

sample_id      <- opt$sample
tumour_name    <- opt$tumour_name
normal_name    <- opt$normal_name
min_depth      <- opt$min_depth
min_alt        <- opt$min_alt_reads
pass_only      <- isTRUE(opt$pass_only)
autosomes_only <- isTRUE(opt$autosomes_only)
max_total_cn   <- opt$max_total_cn
min_snvs       <- opt$min_snvs
min_state_snvs <- opt$min_snvs_per_state
ccf_tol        <- opt$ccf_tolerance
max_exceed     <- opt$max_exceed_fraction
out_prefix     <- if (is.null(opt$out_prefix)) paste0(sample_id, ".vaf_eval") else opt$out_prefix
exceed_p       <- 0.01
peak_rel_height <- 0.25

chrom_levels <- c(as.character(1:22), "X", "Y")
keep_chroms <- if (autosomes_only) as.character(1:22) else c(as.character(1:22), "X")

# Battenberg profile
pp <- read.delim(opt$purity_ploidy)
purity <- as.numeric(pp$purity[1])
ploidy <- as.numeric(pp$ploidy[1])
if (!is.finite(purity) || purity <= 0) stop("Invalid Battenberg purity: ", pp$purity[1])

cn <- read.delim(opt$copynumber, check.names = FALSE)
required_cols <- c("chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A")
if (!all(required_cols %in% colnames(cn))) {
  stop("Copy number file lacks columns: ", paste(setdiff(required_cols, colnames(cn)), collapse = ", "))
}
for (col in c("nMaj2_A", "nMin2_A", "frac2_A")) if (!col %in% colnames(cn)) cn[[col]] <- NA
cn <- cn %>%
  mutate(
    chr      = sub("^chr", "", as.character(chr)),
    startpos = as.numeric(startpos),
    endpos   = as.numeric(endpos),
    nMaj     = as.numeric(nMaj1_A),
    nMin     = as.numeric(nMin1_A),
    clonal   = !is.na(nMaj) & !is.na(nMin) &
               (is.na(suppressWarnings(as.numeric(nMaj2_A))) | is.na(suppressWarnings(as.numeric(frac2_A))) | as.numeric(frac1_A) >= 0.999),
    total_cn = nMaj + nMin,
    state    = ifelse(clonal, paste0(nMaj, "+", nMin), "subclonal"),
    seg_id   = row_number()
  )
seg_len <- cn$endpos - cn$startpos
frac_genome_clonal <- sum(seg_len[cn$clonal], na.rm = TRUE) / sum(seg_len, na.rm = TRUE)

expected_vaf <- function(m, n_total, p) m * p / (p * n_total + 2 * (1 - p))

# VCF parsing
pick_tumour_column <- function(samples) {
  if (tumour_name %in% samples) return(tumour_name)
  for (cand in c("TUMOR", "TUMOUR", "tumor", "tumour")) if (cand %in% samples) return(cand)
  if (length(samples) == 1) return(samples)
  if (normal_name %in% samples && length(samples) == 2) return(setdiff(samples, normal_name))
  hit <- grep(tumour_name, samples, fixed = TRUE, value = TRUE)
  if (nzchar(tumour_name) && length(hit) == 1) return(hit)
  warning("Could not identify the tumour column among: ", paste(samples, collapse = ", "),
          "; using the last sample column")
  samples[length(samples)]
}

empty_snvs <- data.frame(chr = character(), pos = numeric(), ref = character(), alt = character(),
                         ref_count = numeric(), alt_count = numeric())

read_vcf_snvs <- function(path) {
  vcf <- read.vcfR(path, verbose = FALSE)
  stats <- list(n_records = nrow(vcf@fix), tumour_column = NA_character_)
  if (nrow(vcf@fix) == 0 || ncol(vcf@gt) < 2) return(list(snvs = empty_snvs, stats = stats))

  tcol <- pick_tumour_column(colnames(vcf@gt)[-1])
  stats$tumour_column <- tcol

  fix <- getFIX(vcf)
  chr <- sub("^chr", "", fix[, "CHROM"])
  ref <- toupper(fix[, "REF"])
  alt <- toupper(fix[, "ALT"])
  keep <- ref %in% c("A", "C", "G", "T") & alt %in% c("A", "C", "G", "T") & chr %in% keep_chroms
  # vcfR reads a missing FILTER (".") as NA
  if (pass_only) keep <- keep & (is.na(fix[, "FILTER"]) | fix[, "FILTER"] %in% c("PASS", "."))
  stats$n_snvs_pass <- sum(keep)
  if (!any(keep)) return(list(snvs = empty_snvs, stats = stats))

  vcf <- vcf[keep, ]
  v <- data.frame(chr = chr[keep], pos = as.numeric(getPOS(vcf)), ref = ref[keep], alt = alt[keep])
  n <- nrow(v)
  # tumour value of a FORMAT field; NA where the record lacks it
  field <- function(key) extract.gt(vcf, element = key, as.numeric = FALSE, IDtoRowNames = FALSE)[, tcol]
  value_at <- function(x, k) suppressWarnings(as.numeric(vapply(strsplit(x, ",", fixed = TRUE),
    function(z) if (length(z) >= k) z[k] else NA_character_, character(1))))

  fmt <- strsplit(vcf@gt[, "FORMAT"], ":", fixed = TRUE)
  is_strelka <- vapply(fmt, function(f) all(c("AU", "CU", "GU", "TU") %in% f), logical(1))
  has_ad     <- vapply(fmt, function(f) "AD" %in% f, logical(1)) & !is_strelka
  has_af_dp  <- vapply(fmt, function(f) all(c("AF", "DP") %in% f), logical(1)) & !is_strelka & !has_ad

  ref_count <- rep(NA_real_, n)
  alt_count <- rep(NA_real_, n)
  if (any(is_strelka)) {
    # Strelka2 per-base tier 1 counts (AU, CU, GU, TU)
    tier1 <- matrix(vapply(c("A", "C", "G", "T"), function(b) value_at(field(paste0(b, "U")), 1), numeric(n)),
                    nrow = n, dimnames = list(NULL, c("A", "C", "G", "T")))
    i <- which(is_strelka)
    ref_count[i] <- tier1[cbind(i, match(v$ref[i], colnames(tier1)))]
    alt_count[i] <- tier1[cbind(i, match(v$alt[i], colnames(tier1)))]
  }
  if (any(has_ad)) {
    i <- which(has_ad)
    ad <- field("AD")[i]
    ref_count[i] <- value_at(ad, 1)
    alt_count[i] <- value_at(ad, 2)
  }
  if (any(has_af_dp)) {
    i <- which(has_af_dp)
    dp <- value_at(field("DP")[i], 1)
    af <- value_at(field("AF")[i], 1)
    alt_count[i] <- round(af * dp)
    ref_count[i] <- dp - alt_count[i]
  }
  stats$count_source <- paste(c(if (any(is_strelka)) "strelka_tier1", if (any(has_ad)) "AD",
                                if (any(has_af_dp)) "AF*DP"), collapse = ",")
  list(snvs = data.frame(chr = v$chr, pos = v$pos, ref = v$ref, alt = v$alt,
                         ref_count = ref_count, alt_count = alt_count),
       stats = stats)
}

vcf <- read_vcf_snvs(opt$vcf)
snv <- vcf$snvs %>%
  mutate(depth = ref_count + alt_count, vaf = alt_count / depth) %>%
  filter(is.finite(vaf), depth >= min_depth, alt_count >= min_alt)
cat(sprintf("[evaluate_vaf] %s: %s VCF records, %s PASS SNVs, %d SNVs after depth filters (tumour column: %s, counts: %s)\n",
            sample_id, vcf$stats$n_records, ifelse(is.null(vcf$stats$n_snvs_pass), 0, vcf$stats$n_snvs_pass),
            nrow(snv), vcf$stats$tumour_column, ifelse(is.null(vcf$stats$count_source), "NA", vcf$stats$count_source)))

# Assign SNVs to segments
snv$seg_id <- rep(NA_integer_, nrow(snv))
for (ch in unique(snv$chr)) {
  s <- which(cn$chr == ch)
  if (length(s) == 0) next
  s <- s[order(cn$startpos[s])]
  i <- which(snv$chr == ch)
  j <- findInterval(snv$pos[i], cn$startpos[s])
  ok <- j > 0
  cand <- s[j[ok]]
  inside <- snv$pos[i][ok] <= cn$endpos[cand]
  snv$seg_id[i[ok][inside]] <- cand[inside]
}

snv <- snv %>%
  left_join(cn %>% select(seg_id, clonal, nMaj, nMin, total_cn, state), by = "seg_id") %>%
  mutate(
    on_clonal_segment = !is.na(seg_id) & clonal %in% TRUE,
    assessed          = on_clonal_segment & nMaj >= 1 & total_cn <= max_total_cn,
    denom             = purity * total_cn + 2 * (1 - purity),
    expected_vaf_m1   = ifelse(assessed, purity / denom, NA_real_),
    expected_vaf_max  = ifelse(assessed, nMaj * purity / denom, NA_real_),
    multiplicity      = ifelse(assessed, vaf * denom / purity, NA_real_),
    best_m            = ifelse(assessed, pmin(pmax(round(multiplicity), 1), nMaj), NA_real_),
    ccf               = multiplicity / best_m,
    p_exceed          = ifelse(assessed,
                               pbinom(alt_count - 1, depth, pmin(expected_vaf_max, 1), lower.tail = FALSE),
                               NA_real_),
    exceeds_max       = p_exceed < exceed_p
  )
clonal_snv <- snv %>% filter(assessed)

# Peak finding
kde_peaks <- function(x, from, to) {
  x <- x[is.finite(x) & x >= from & x <= to]
  if (length(x) < 10) return(data.frame(position = numeric(0), height = numeric(0)))
  d <- density(x, from = from, to = to, n = 512)
  y <- d$y
  idx <- which(diff(sign(diff(y))) == -2) + 1
  if (y[length(y)] > y[length(y) - 1]) idx <- c(idx, length(y))   # peak at the upper bound
  pk <- data.frame(position = d$x[idx], height = y[idx])
  pk[pk$height >= peak_rel_height * max(y), , drop = FALSE]
}
weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

# Per clonal copy number state
state_rows <- list()
for (st in unique(clonal_snv$state)) {
  d <- clonal_snv[clonal_snv$state == st, ]
  n_maj <- d$nMaj[1]; n_min <- d$nMin[1]; n_tot <- n_maj + n_min
  exp_peaks <- expected_vaf(seq_len(n_maj), n_tot, purity)
  peak_vaf <- NA_real_; nearest_m <- NA_real_; implied <- NA_real_; obs_peaks <- ""
  if (nrow(d) >= min_state_snvs) {
    pk <- kde_peaks(d$vaf, 0.02, 1)
    if (nrow(pk) > 0) {
      obs_peaks <- paste(sprintf("%.3f", pk$position), collapse = ";")
      peak_vaf <- max(pk$position)                                   # rightmost prominent peak
      nearest_m <- which.min(abs(exp_peaks - peak_vaf))
      # purity implied by the peak: solve VAF = m p / (p N + 2 (1 - p)) for p
      denom_p <- nearest_m - peak_vaf * (n_tot - 2)
      implied <- if (denom_p > 0) 2 * peak_vaf / denom_p else NA_real_
    }
  }
  state_rows[[st]] <- data.frame(
    sample = sample_id, state = st, nMaj = n_maj, nMin = n_min, total_cn = n_tot,
    n_snvs = nrow(d),
    expected_vaf_peaks = paste(sprintf("%.3f", exp_peaks), collapse = ";"),
    observed_vaf_peaks = obs_peaks,
    clonal_peak_vaf = round(peak_vaf, 4),
    nearest_multiplicity = nearest_m,
    expected_vaf_nearest = if (is.na(nearest_m)) NA_real_ else round(exp_peaks[nearest_m], 4),
    implied_purity = round(implied, 4),
    median_ccf = round(median(d$ccf, na.rm = TRUE), 4),
    frac_exceed_max = round(mean(d$exceeds_max, na.rm = TRUE), 4)
  )
}
states <- if (length(state_rows)) bind_rows(state_rows) %>% arrange(desc(n_snvs)) else
  data.frame(sample = character(), state = character(), nMaj = numeric(), nMin = numeric(), total_cn = numeric(),
             n_snvs = integer(), expected_vaf_peaks = character(), observed_vaf_peaks = character(),
             clonal_peak_vaf = numeric(), nearest_multiplicity = numeric(), expected_vaf_nearest = numeric(),
             implied_purity = numeric(), median_ccf = numeric(), frac_exceed_max = numeric())

# Sample-level metrics
n_clonal <- nrow(clonal_snv)
ccf_peaks <- kde_peaks(clonal_snv$ccf, 0.2, 3)
ccf_peak <- if (nrow(ccf_peaks) > 0) max(ccf_peaks$position) else NA_real_
frac_ccf_near_1 <- if (n_clonal > 0) mean(abs(clonal_snv$ccf - 1) <= 0.25, na.rm = TRUE) else NA_real_
frac_exceed <- if (n_clonal > 0) mean(clonal_snv$exceeds_max, na.rm = TRUE) else NA_real_
snv_purity <- weighted_median(states$implied_purity, states$n_snvs)
multi <- clonal_snv %>% filter(nMaj >= 2)
frac_m_gt1 <- if (nrow(multi) > 0) mean(multi$best_m >= 2) else NA_real_

reasons <- character(0)
if (n_clonal < min_snvs) {
  qc_flag <- "INSUFFICIENT_SNVS"
  reasons <- sprintf("%d SNVs on clonal segments (< %d)", n_clonal, as.integer(min_snvs))
} else {
  if (is.na(ccf_peak) || abs(ccf_peak - 1) > ccf_tol) {
    reasons <- c(reasons, sprintf("clonal CCF peak %.2f outside 1 +/- %.2f", ccf_peak, ccf_tol))
  }
  if (is.finite(frac_exceed) && frac_exceed > max_exceed) {
    reasons <- c(reasons, sprintf("%.1f%% of SNVs above maximum expected VAF", 100 * frac_exceed))
  }
  qc_flag <- if (length(reasons) == 0) "PASS" else "REVIEW"
}

metrics <- data.frame(
  sample                    = sample_id,
  tumour_name               = tumour_name,
  purity                    = purity,
  ploidy                    = ploidy,
  frac_genome_clonal        = round(frac_genome_clonal, 4),
  n_vcf_records             = vcf$stats$n_records,
  n_snvs_filtered           = nrow(snv),
  n_snvs_clonal_segments    = n_clonal,
  n_snvs_subclonal_segments = sum(!is.na(snv$seg_id) & !(snv$clonal %in% TRUE)),
  n_snvs_unassigned         = sum(is.na(snv$seg_id)),
  n_states_evaluated        = sum(states$n_snvs >= min_state_snvs),
  ccf_clonal_peak           = round(ccf_peak, 4),
  frac_ccf_near_1           = round(frac_ccf_near_1, 4),
  frac_exceed_max           = round(frac_exceed, 4),
  snv_implied_purity        = round(snv_purity, 4),
  purity_delta              = round(snv_purity - purity, 4),
  frac_multiplicity_gt1     = round(frac_m_gt1, 4),
  qc_flag                   = qc_flag,
  qc_reason                 = if (length(reasons)) paste(reasons, collapse = "; ") else "",
  vcf_tumour_column         = vcf$stats$tumour_column
)

write.table(metrics, paste0(out_prefix, ".metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(states, paste0(out_prefix, ".states.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
gz <- gzfile(paste0(out_prefix, ".snvs.tsv.gz"), "w")
write.table(snv %>% mutate(sample = sample_id, .before = 1) %>%
              select(sample, chr, pos, ref, alt, ref_count, alt_count, depth, vaf, seg_id, state, nMaj, nMin,
                     on_clonal_segment, expected_vaf_m1, expected_vaf_max, multiplicity, best_m, ccf, p_exceed),
            gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
close(gz)

# Plots
title_txt <- sprintf("%s  |  Battenberg purity %.2f, ploidy %.2f  |  VAF QC: %s", sample_id, purity, ploidy, qc_flag)
pdf(paste0(out_prefix, ".pdf"), width = 12, height = 8)

if (n_clonal == 0) {
  plot.new()
  title(main = title_txt)
  text(0.5, 0.5, "No SNVs passed filters on clonal copy number segments")
} else {
  # VAF per clonal copy number state with expected clonal peaks
  top <- states %>% filter(n_snvs >= 10) %>% head(12)
  if (nrow(top) > 0) {
    lab <- setNames(sprintf("%s  (N=%d, n=%d)", top$state, as.integer(top$total_cn), as.integer(top$n_snvs)), top$state)
    d1 <- clonal_snv %>% filter(state %in% top$state) %>% mutate(panel = factor(lab[state], levels = lab))
    expl <- bind_rows(lapply(seq_len(nrow(top)), function(k) {
      data.frame(panel = lab[[top$state[k]]], m = seq_len(top$nMaj[k]),
                 vaf = expected_vaf(seq_len(top$nMaj[k]), top$total_cn[k], purity))
    })) %>% mutate(panel = factor(panel, levels = lab), m = factor(m))
    p1 <- ggplot(d1, aes(x = vaf)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 0.02, boundary = 0, fill = "grey75", colour = NA) +
      geom_density(colour = "grey20") +
      geom_vline(data = expl, aes(xintercept = vaf, colour = m), linetype = "dashed") +
      facet_wrap(~ panel, scales = "free_y") +
      coord_cartesian(xlim = c(0, 1)) +
      labs(title = title_txt,
           subtitle = "VAF of SNVs on clonal segments by copy number state (nMaj+nMin); dashed lines: expected clonal VAF for multiplicity m",
           x = "Variant allele frequency", y = "Density", colour = "m") +
      theme_bw(base_size = 10)
    print(p1)
  }

  # Pooled CCF
  p2 <- ggplot(clonal_snv %>% filter(ccf <= 3), aes(x = ccf)) +
    geom_histogram(binwidth = 0.05, boundary = 0, fill = "steelblue", colour = NA) +
    geom_vline(xintercept = 1, colour = "black") +
    annotate("rect", xmin = 1 - ccf_tol, xmax = 1 + ccf_tol, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "darkgreen") +
    { if (is.finite(ccf_peak)) geom_vline(xintercept = ccf_peak, colour = "red", linetype = "dashed") } +
    labs(title = title_txt,
         subtitle = sprintf("Cancer cell fraction (best multiplicity) of %d SNVs on clonal segments; red: clonal peak %.2f; shaded: 1 +/- %.2f",
                            n_clonal, ccf_peak, ccf_tol),
         x = "CCF", y = "SNVs") +
    theme_bw(base_size = 10)
  print(p2)

  # Genome-wide VAF with expected clonal VAF per segment
  segs <- cn %>%
    filter(clonal, nMaj >= 1, chr %in% keep_chroms) %>%
    mutate(exp_m1 = expected_vaf(1, total_cn, purity), exp_max = expected_vaf(nMaj, total_cn, purity),
           chr = factor(chr, levels = chrom_levels))
  d3 <- snv %>% filter(!is.na(seg_id)) %>%
    mutate(chr = factor(chr, levels = chrom_levels),
           segment = ifelse(on_clonal_segment, "clonal segment", "subclonal segment"))
  p3 <- ggplot() +
    geom_point(data = d3, aes(x = pos / 1e6, y = vaf, colour = segment), size = 0.3, alpha = 0.4) +
    geom_segment(data = segs, aes(x = startpos / 1e6, xend = endpos / 1e6, y = exp_m1, yend = exp_m1),
                 colour = "red", linewidth = 0.6) +
    geom_segment(data = segs %>% filter(nMaj >= 2), aes(x = startpos / 1e6, xend = endpos / 1e6, y = exp_max, yend = exp_max),
                 colour = "darkorange", linewidth = 0.6) +
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = c("clonal segment" = "grey30", "subclonal segment" = "purple")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title_txt, subtitle = "Expected clonal VAF: red m = 1, orange m = nMaj",
         x = "Position (Mb)", y = "VAF", colour = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing = unit(0.05, "lines"),
          legend.position = "bottom")
  print(p3)
}
invisible(dev.off())

print(metrics)
