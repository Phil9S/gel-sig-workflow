#!/usr/bin/env Rscript
# Paired tumour/normal WGS Battenberg run, based on inst/example/battenberg_wgs.R.
# Works with Battenberg v3.x and v2.2.x; outputs are normalised to the v3 names
# (<samplename>_copynumber.txt, _copynumber_extended.txt, _purity_ploidy.txt).

suppressPackageStartupMessages(library(Battenberg))
options(bitmapType = "cairo", stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    if (!startsWith(args[i], "--")) stop("Unexpected argument: ", args[i])
    key <- sub("^--", "", args[i])
    if (i == length(args) || startsWith(args[i + 1], "--")) {
      out[[key]] <- "TRUE"
      i <- i + 1
    } else {
      out[[key]] <- args[i + 1]
      i <- i + 2
    }
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
req <- function(key) {
  if (is.null(opt[[key]]) || opt[[key]] == "") stop("Missing required argument --", key)
  opt[[key]]
}
num <- function(key, default = NULL) if (is.null(opt[[key]])) default else as.numeric(opt[[key]])
flag <- function(key, default = FALSE) if (is.null(opt[[key]])) default else toupper(opt[[key]]) %in% c("TRUE", "T", "1", "YES")

ref_dir <- req("ref_dir")
ref_path <- function(key) {
  x <- opt[[key]]
  if (is.null(x) || x %in% c("", "null", "NA")) return(NULL)
  if (startsWith(x, "/")) x else file.path(ref_dir, x)
}

samplename <- req("samplename")
normalname <- req("normalname")
sex <- toupper(req("sex"))
is_male <- sex %in% c("XY", "MALE", "M")
genomebuild <- req("genomebuild")

cat("[run_battenberg] Battenberg version:", as.character(packageVersion("Battenberg")), "\n")
cat("[run_battenberg] sample:", samplename, " normal:", normalname, " sex:", sex, " build:", genomebuild, "\n")

# impute_info.txt
# Rewrite reference paths to resolve inside ref_dir. Missing absolute paths (bundle
# prepared elsewhere) are re-anchored by matching their trailing path components.
relocate <- function(p) {
  if (is.na(p) || p == "" || file.exists(p)) return(p)
  parts <- strsplit(p, "/", fixed = TRUE)[[1]]
  parts <- parts[parts != ""]
  for (k in seq_len(min(length(parts), 4))) {
    cand <- file.path(ref_dir, paste(tail(parts, k), collapse = "/"))
    if (file.exists(cand)) return(cand)
  }
  if (!startsWith(p, "/")) return(file.path(ref_dir, p))
  p
}
impute_info_in <- ref_path("impute_info")
if (is.null(impute_info_in) || !file.exists(impute_info_in)) stop("impute_info file not found: ", impute_info_in)
# whitespace separated, as read by Battenberg::parse.imputeinfofile
impute_info <- read.table(impute_info_in, header = FALSE, stringsAsFactors = FALSE, quote = "")
for (col in 2:4) impute_info[[col]] <- vapply(impute_info[[col]], relocate, character(1), USE.NAMES = FALSE)
missing_ref <- unique(unlist(impute_info[, 2:4]))
missing_ref <- missing_ref[!file.exists(missing_ref)]
if (length(missing_ref) > 0) {
  warning("impute_info references files that could not be found, e.g.: ", paste(head(missing_ref, 3), collapse = ", "))
}
impute_info_local <- "impute_info.local.txt"
write.table(impute_info, impute_info_local, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# battenberg() call
args <- list(
  analysis                    = "paired",
  samplename                  = samplename,
  normalname                  = normalname,
  sample_data_file            = req("tumour_bam"),
  normal_data_file            = req("normal_bam"),
  ismale                      = is_male,
  imputeinfofile              = impute_info_local,
  g1000prefix                 = ref_path("g1000_loci_prefix"),
  g1000allelesprefix          = ref_path("g1000_alleles_prefix"),
  gccorrectprefix             = ref_path("gc_prefix"),
  repliccorrectprefix         = ref_path("rt_prefix"),
  problemloci                 = ref_path("problem_loci"),
  chrom_coord_file            = ref_path("chrom_coord_file"),
  data_type                   = "wgs",
  impute_exe                  = "impute2",
  allelecounter_exe           = "alleleCounter",
  javajre                     = "java",
  usebeagle                   = flag("use_beagle", TRUE),
  beaglejar                   = ref_path("beagle_jar"),
  beagleref.template          = ref_path("beagle_ref_template"),
  beagleplink.template        = ref_path("beagle_plink_template"),
  beaglemaxmem                = num("beagle_max_mem", 10),
  beaglenthreads              = 1,
  beaglewindow                = num("beagle_window", 40),
  beagleoverlap               = num("beagle_overlap", 4),
  nthreads                    = num("cpus", 8),
  platform_gamma              = num("platform_gamma", 1),
  phasing_gamma               = num("phasing_gamma", 1),
  segmentation_gamma          = num("segmentation_gamma", 10),
  segmentation_kmin           = num("segmentation_kmin", 3),
  phasing_kmin                = num("phasing_kmin", 1),
  clonality_dist_metric       = num("clonality_dist_metric", 0),
  ascat_dist_metric           = num("ascat_dist_metric", 1),
  min_ploidy                  = num("min_ploidy", 1.6),
  max_ploidy                  = num("max_ploidy", 4.8),
  min_rho                     = num("min_rho", 0.1),
  max_rho                     = num("max_rho", 1.0),
  min_goodness                = num("min_goodness", 0.63),
  uninformative_BAF_threshold = num("uninformative_baf_threshold", 0.51),
  min_normal_depth            = num("min_normal_depth", 10),
  min_base_qual               = num("min_base_qual", 20),
  min_map_qual                = num("min_map_qual", 35),
  calc_seg_baf_option         = num("calc_seg_baf_option", 3),
  max_allowed_state           = num("max_allowed_state", 250),
  skip_allele_counting        = FALSE,
  skip_preprocessing          = FALSE,
  skip_phasing                = FALSE,
  genomebuild                 = genomebuild,
  enhanced_grid_search        = flag("enhanced_grid_search", FALSE)
)
args <- args[!vapply(args, is.null, logical(1))]

# chrom_coord_file is only used by the cell_line/germline modes; drop it when absent
if (!is.null(args$chrom_coord_file) && !file.exists(args$chrom_coord_file)) {
  cat("[run_battenberg] chrom_coord_file not found, not used for paired analysis:", args$chrom_coord_file, "\n")
  args$chrom_coord_file <- NULL
}
for (key in c("g1000prefix", "g1000allelesprefix", "gccorrectprefix", "repliccorrectprefix")) {
  if (!is.null(args[[key]]) && length(Sys.glob(paste0(args[[key]], "*"))) == 0) {
    stop("No reference files match ", key, " = ", args[[key]], "*  (check params.genomes in conf/genomes.config)")
  }
}
if (!is.null(args$problemloci) && !file.exists(args$problemloci)) stop("problem loci file not found: ", args$problemloci)
if (isTRUE(args$usebeagle)) {
  if (is.null(args$beaglejar) || !file.exists(args$beaglejar)) stop("Beagle jar not found: ", args$beaglejar)
  for (key in c("beagleref.template", "beagleplink.template")) {
    chr1 <- sub("CHROMNAME", "1", args[[key]], fixed = TRUE)
    if (is.null(args[[key]]) || !file.exists(chr1)) stop("Beagle reference not found for chromosome 1: ", chr1, " (from ", key, ")")
  }
}

# Map onto the installed signature (v2.2.x names differ from v3.x)
bb_formals <- names(formals(Battenberg::battenberg))
renames <- list(
  samplename           = c("tumourname"),
  sample_data_file     = c("tumour_data_file"),
  beagleref.template   = c("beagleref"),
  beagleplink.template = c("beagleplink"),
  genomebuild          = c("GENOMEBUILD")
)
for (new_name in names(renames)) {
  if (new_name %in% names(args) && !(new_name %in% bb_formals)) {
    alt <- intersect(renames[[new_name]], bb_formals)
    if (length(alt) > 0) names(args)[names(args) == new_name] <- alt[1]
  }
}
unsupported <- setdiff(names(args), bb_formals)
if (length(unsupported) > 0) {
  cat("[run_battenberg] not supported by this Battenberg version, ignored:", paste(unsupported, collapse = ", "), "\n")
  args <- args[names(args) %in% bb_formals]
}

cat("[run_battenberg] battenberg() arguments:\n")
str(args)
do.call(Battenberg::battenberg, args)

# Normalise outputs to v3 names
cn_file <- paste0(samplename, "_copynumber.txt")
cn_ext_file <- paste0(samplename, "_copynumber_extended.txt")
pp_file <- paste0(samplename, "_purity_ploidy.txt")
rp_file <- paste0(samplename, "_rho_and_psi.txt")
cn_cols <- c("chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")

if (!file.exists(cn_file)) {
  old <- paste0(samplename, "_subclones.txt")
  if (!file.exists(old)) stop("Battenberg finished without writing ", cn_file, " or ", old)
  cat("[run_battenberg] converting", old, "to", cn_file, "\n")
  sc <- read.delim(old, check.names = FALSE)
  write.table(sc, cn_ext_file, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sc[, intersect(cn_cols, colnames(sc))], cn_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
if (!file.exists(cn_ext_file)) file.copy(cn_file, cn_ext_file)

rho_psi <- if (file.exists(rp_file)) read.table(rp_file, header = TRUE, sep = "\t") else NULL
frac_genome <- function(col) {
  if (is.null(rho_psi) || !("FRAC_GENOME" %in% rownames(rho_psi)) || !(col %in% colnames(rho_psi))) return(NA_real_)
  as.numeric(rho_psi[["FRAC_GENOME", col]])
}

pp <- NULL
if (file.exists(pp_file)) {
  pp <- read.delim(pp_file)
} else if (file.exists(paste0(samplename, "_cellularity_ploidy.txt"))) {
  pp <- read.delim(paste0(samplename, "_cellularity_ploidy.txt"))
  colnames(pp)[colnames(pp) == "cellularity"] <- "purity"
}
if (is.null(pp)) pp <- data.frame(purity = frac_genome("rho"), ploidy = frac_genome("ploidy"))
if (!("psi" %in% colnames(pp))) pp$psi <- frac_genome("psi")
write.table(pp[, c("purity", "ploidy", "psi")], pp_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("[run_battenberg] purity/ploidy:\n")
print(pp)
sessionInfo()
