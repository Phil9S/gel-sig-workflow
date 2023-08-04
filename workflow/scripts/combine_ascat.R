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

## compute total CN
segment.table$segVal <- segment.table$nAraw + segment.table$nBraw
segment.table <- segment.table[,c("chr","startpos","endpos","segVal","sample")]
colnames(segment.table) <- c("chromosome","start","end","segVal","sample")

## Set smoothing factor
SMOOTHING_FACTOR <- 0.12

## Smooth segments within SMOOTHING_FACTOR of adjacent segments
segment.table.smoothed.rounded <- segment.table %>%
  #mutate(segVal = round(segVal,digits = 2)) %>%
  group_by(chromosome,sample) %>%
  mutate(seg_diff = abs(segVal - lag(segVal))) %>%
  mutate(chng = ifelse(seg_diff > SMOOTHING_FACTOR,"TRUE","FALSE")) %>%
  mutate(chng = as.logical(ifelse(is.na(chng),"TRUE",chng))) %>%
  mutate(comb = cumsum(chng)) %>%
  group_by(chromosome,sample,comb) %>%
  select(-chng) %>%
  summarise(across(start,min),across(end,max),across(segVal,median)) %>%
  select(chromosome,start,end,segVal,sample) %>%
  mutate(chromosome = factor(chromosome,levels=c(1:22,"X","Y"))) %>%
  arrange(sample,chromosome,start)

#print(head(segment.table.smoothed.rounded))

## Export segment table
write.table(x = segment.table.smoothed.rounded,
	file = OUTPUT,
	quote = F,
	row.names = F,
	col.names = T,
	sep = "\t")

