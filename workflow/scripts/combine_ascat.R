#In singularity
## allow cmdline args
args = commandArgs(trailingOnly=TRUE)

## Load libraries
library(dplyr)

## Assign cmdlind args
INPUT <- args[2:length(args)]
OUTPUT <- args[1]

## Load seg tables
file.list <- lapply(INPUT,FUN=function(x){read.table(x,sep="\t",header=TRUE)})

## Compute total and collapse seg tables
segment.table <- do.call(rbind,lapply(file.list,FUN=function(x){
	return(as.data.frame(x))
}))

## remove failed arrays
segment.table <- segment.table[!is.na(segment.table$chr),]

## Export segment table
write.table(x = segment.table,
	file = OUTPUT,
	quote = F,
	row.names = F,
	col.names = T,
	sep = "\t")

