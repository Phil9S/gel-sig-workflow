# In singularity

library(ASCAT)

p <- sessionInfo()

print(snakemake@output)

writeLines(p,snakemake@output)
