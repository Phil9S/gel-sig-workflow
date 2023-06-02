# In singularity

args = commandArgs(trailingOnly=TRUE)

INPUT <- args[1]
OUTPUT <- args[2]
GENOME <- args[3]
METHOD <- args[4]
THREADS <- args[5]

library(CINSignatureQuantification)

experiment <- paste0("quantifyCNSignatures-",format(Sys.time(), '%Y-%m-%d %H:%M:%S'))

quant.sigs <- quantifyCNSignatures(object=INPUT,
    experimentName=experiment,
    method=METHOD,
    cores=THREADS,
    build=GENOME)

saveRDS(object=quant.sigs,file=OUTPUT)
