# In singularity

args = commandArgs(trailingOnly=TRUE)

library(ASCAT)

# Sample args
SAMPLE <- args[1]
TUMOUR <- args[2]
NORMAL <- args[3]
TUMOUR_NAME <- args[4]
NORMAL_NAME <- args[5]
SEX <- args[6]
CANCER <- args[7]

# Build args
GENOME <- args[8]

# out args
OUTPUT <- args[9]
OUT_DIR <- args[10]
OUTTARGET <- paste0(OUT_DIR,"fitting/",CANCER,"/",SAMPLE,"/")
OUTTARGET_T <- paste0(OUT_DIR,"fitting/",CANCER,"/",SAMPLE,"/",TUMOUR_NAME)
OUTTARGET_N <- paste0(OUT_DIR,"fitting/",CANCER,"/",SAMPLE,"/",NORMAL_NAME)

# resource args
THREADS <- args[11]

# Sing args
ALLELECOUNTPATH <- "/usr/local/bin/alleleCounter"
ALLELESPREFIX <- paste0("/ref/",GENOME,"/G1000_alleles_",GENOME,"_chr")
LOCIPREFIX <- paste0("/ref/",GENOME,"/G1000_loci_",GENOME,"_chr")
GCFILE <- paste0("/ref/",GENOME,"/GC_G1000_",GENOME,".txt")
RTFILE <- paste0("/ref/",GENOME,"/RT_G1000_",GENOME,".txt")

#print(ALLELESPREFIX)
#print(LOCIPREFIX)

# Manually adjust the prepareHTS function to set these to prevent output from alleleCounter dumping in $PWD
## tumourAlleleCountsFile.prefix 
## normalAlleleCountsFile.prefix
ALLELECOUNTSPREFIX <- paste0(OUTTARGET,"ALLELE_COUNTS/")

# Make dir if not present
if(!dir.exists(ALLELECOUNTSPREFIX)){
  dir.create(ALLELECOUNTSPREFIX)
}
# Source modified function to use the adjusted output locs
source("workflow/scripts/ascat.prepareHTS_modified.R")

ascat.prepareHTS(
  tumourseqfile = TUMOUR,
  normalseqfile = NORMAL,
  tumourname = TUMOUR_NAME,
  normalname = NORMAL_NAME,
  allelecounter_exe = ALLELECOUNTPATH,
  alleles.prefix = ALLELESPREFIX,
  loci.prefix = LOCIPREFIX,
  gender = SEX,
  genomeVersion = GENOME,
  nthreads = THREADS,
  tumourLogR_file = paste0(OUTTARGET_T,"_tumor_LogR.txt"),
  tumourBAF_file = paste0(OUTTARGET_T,"_tumor_BAF.txt"),
  normalLogR_file = paste0(OUTTARGET_N,"_germline_LogR.txt"),
  normalBAF_file = paste0(OUTTARGET_N,"_germline_BAF.txt"))

ascat.bc = ascat.loadData(
    Tumor_LogR_file = paste0(OUTTARGET_T,"_tumor_LogR.txt"),
    Tumor_BAF_file = paste0(OUTTARGET_T,"_tumor_BAF.txt"),
    Germline_LogR_file = paste0(OUTTARGET_N,"_germline_LogR.txt"),
    Germline_BAF_file = paste0(OUTTARGET_N,"_germline_BAF.txt"),
    gender = SEX,
    genomeVersion = GENOME)

ascat.plotRawData(ascat.bc,img.dir=OUTTARGET, img.prefix = "Before_correction_")

ascat.bc = ascat.correctLogR(ascat.bc,
    GCcontentfile = GCFILE,
    replictimingfile = RTFILE)

ascat.plotRawData(ascat.bc,img.dir=OUTTARGET, img.prefix = "After_correction_")

ascat.bc = ascat.aspcf(ascat.bc,out.dir=OUTTARGET)

ascat.plotSegmentedData(ascat.bc,img.dir=OUTTARGET)

ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T,img.dir=OUTTARGET)

QC = ascat.metrics(ascat.bc,ascat.output)

save(ascat.bc, ascat.output, QC, file = paste0(OUTTARGET,SAMPLE,"_ASCAT_objects.Rdata"))
write.table(x = QC, file = paste0(OUTTARGET,SAMPLE,"_QC_metrics.tsv"),append = FALSE,
	quote = FALSE, row.names = FALSE, col.names = TRUE, sep ="\t")


if(length(ascat.output$failedarrays) > 0){
	d <- data.frame(sample=c(ascat.output$failedarrays),
		chr=c(NA),startpos=c(NA),endpos=c(NA),
		nMajor=c(NA),nMinor=c(NA),nAraw=c(NA),nBraw=c(NA))

	write.table(x = d,file = paste0(OUTTARGET_T,".segments_raw.txt"),append = FALSE,
			quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
} else {

	sample.adjust <- read.table(paste0(OUTTARGET_T,".segments_raw.txt"),header=TRUE,sep="\t")
	sample.adjust$sample <- rep(SAMPLE,times=nrow(sample.adjust))
	write.table(sample.adjust,file=paste0(OUTTARGET_T,".segments_raw.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
}

sessionInfo()

