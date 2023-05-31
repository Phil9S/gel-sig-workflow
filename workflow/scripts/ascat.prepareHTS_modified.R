# This function is modified to customise the output of allelecounter
ascat.prepareHTS = function(tumourseqfile, normalseqfile, tumourname, normalname, allelecounter_exe, alleles.prefix, loci.prefix, gender, genomeVersion,
                            nthreads=1, tumourLogR_file=NA, tumourBAF_file=NA, normalLogR_file=NA, normalBAF_file=NA, minCounts=10, BED_file=NA,
                            probloci_file=NA, chrom_names=c(1:22,'X'), min_base_qual=20, min_map_qual=35, ref.fasta=NA,
                            skip_allele_counting_tumour=F, skip_allele_counting_normal=F, seed=as.integer(Sys.time())) {
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  doParallel::registerDoParallel(cores=nthreads)
  
  if (is.na(tumourLogR_file)) tumourLogR_file=paste0(tumourname,"_tumourLogR.txt")
  if (is.na(tumourBAF_file)) tumourBAF_file=paste0(tumourname,"_tumourBAF.txt")
  if (is.na(normalLogR_file)) normalLogR_file=paste0(tumourname,"_normalLogR.txt")
  if (is.na(normalBAF_file)) normalBAF_file=paste0(tumourname,"_normalBAF.txt")
  
  if (!skip_allele_counting_tumour) {
    # Obtain allele counts at specific loci for tumour
    foreach::foreach(CHR=chrom_names) %dopar% {
      ascat.getAlleleCounts(seq.file=tumourseqfile,
                            output.file=paste0(ALLELECOUNTSPREFIX,tumourname,"_alleleFrequencies_chr", CHR, ".txt"),
                            loci.file=paste0(loci.prefix, CHR, ".txt"),
                            min.base.qual=min_base_qual,
                            min.map.qual=min_map_qual,
                            allelecounter.exe=allelecounter_exe,
                            ref.fasta=ref.fasta)
    }
  }
  if (!skip_allele_counting_normal) {
    # Obtain allele counts at specific loci for normal
    foreach::foreach(CHR=chrom_names) %dopar% {
      ascat.getAlleleCounts(seq.file=normalseqfile,
                            output.file=paste0(ALLELECOUNTSPREFIX,normalname,"_alleleFrequencies_chr", CHR, ".txt"),
                            loci.file=paste0(loci.prefix, CHR, ".txt"),
                            min.base.qual=min_base_qual,
                            min.map.qual=min_map_qual,
                            allelecounter.exe=allelecounter_exe,
                            ref.fasta=ref.fasta)
    }
  }
  # Obtain BAF and LogR from the raw allele counts
  ascat.getBAFsAndLogRs(samplename=tumourname,
                        tumourAlleleCountsFile.prefix=paste0(ALLELECOUNTSPREFIX,tumourname,"_alleleFrequencies_chr"),
                        normalAlleleCountsFile.prefix=paste0(ALLELECOUNTSPREFIX,normalname,"_alleleFrequencies_chr"),
                        tumourLogR_file=tumourLogR_file,
                        tumourBAF_file=tumourBAF_file,
                        normalLogR_file=normalLogR_file,
                        normalBAF_file=normalBAF_file,
                        alleles.prefix=alleles.prefix,
                        gender=gender,
                        genomeVersion=genomeVersion,
                        chrom_names=chrom_names,
                        minCounts=minCounts,
                        BED_file=BED_file,
                        probloci_file=probloci_file,
                        seed=seed)

  # Synchronise all information
  ascat.synchroniseFiles(samplename=tumourname,
                         tumourLogR_file=tumourLogR_file,
                         tumourBAF_file=tumourBAF_file,
                         normalLogR_file=normalLogR_file,
                         normalBAF_file=normalBAF_file)
}

#' Function to concatenate allele counter output
#' @noRd
readAlleleCountFiles=function(prefix,suffix,chrom_names,minCounts,keep_chr_string=F) {
  files=paste0(prefix,chrom_names,suffix)
  files=files[sapply(files,function(x) file.exists(x) && file.info(x)$size>0)]
  stopifnot(length(files)>0)
  data=do.call(rbind,lapply(files,function(x) {
    tmp=data.frame(data.table::fread(x,sep='\t',showProgress=F,header=T),stringsAsFactors=F)
    tmp=tmp[tmp[,7]>=minCounts,]
    if (nrow(tmp)>0) {
      if (!keep_chr_string) tmp[,1]=gsub('^chr','',tmp[,1])
      rownames(tmp)=paste0(tmp[,1],'_',tmp[,2])
    }
    return(tmp)
  }))
  stopifnot(nrow(data)>0)
  return(data)
}

#' Function to concatenate all alleles files
#' @noRd
readAllelesFiles=function(prefix,suffix,chrom_names,add_chr_string=F) {
  files=paste0(prefix,chrom_names,suffix)
  files=files[sapply(files,function(x) file.exists(x) && file.info(x)$size>0)]
  stopifnot(length(files)>0)
  data=do.call(rbind,lapply(files,function(x) {
    tmp=data.frame(data.table::fread(x,sep='\t',showProgress=F,header=T))
    tmp=tmp[!is.na(tmp[,2] & !is.na(tmp[,3])),]
    tmp=tmp[!duplicated(tmp[,1]),]
    tmp$chromosome=gsub(paste0(prefix,'(',paste(chrom_names,collapse='|'),')',suffix),'\\1',x)
    if (add_chr_string) tmp$chromosome=paste0('chr',tmp$chromosome)
    tmp=tmp[,c(4,1:3)]
    rownames(tmp)=paste0(tmp[,1],'_',tmp[,2])
    return(tmp)
  }))
  stopifnot(nrow(data)>0)
  return(data)
}

#' Function to concatenate all loci files
#' @noRd
readLociFiles=function(prefix,suffix,chrom_names,keep_chr_string=F) {
  files=paste0(prefix,chrom_names,suffix)
  files=files[sapply(files,function(x) file.exists(x) && file.info(x)$size>0)]
  stopifnot(length(files)>0)
  data=do.call(rbind,lapply(files,function(x) {
    tmp=data.frame(data.table::fread(x,sep='\t',showProgress=F,header=F))
    if (!keep_chr_string) tmp[,1]=gsub('^chr','',tmp[,1])
    rownames(tmp)=paste0(tmp[,1],'_',tmp[,2])
    return(tmp)
  }))
  stopifnot(nrow(data)>0)
  return(data)
}
