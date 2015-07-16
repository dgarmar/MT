
##
## Index genotype file (if not done externally before)
##

## WD
setwd("/nfs/users/rg/dgarrido/run-GTEx/bin")

## Look at the paths in the specified order when loading a library
.libPaths(c("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1",
            "/nfs/software/R/packages"))

## Command line input
args <- commandArgs(TRUE)

library(sQTLseekeR)

## Input file: genotype (required format)
genotype.f = args[1]

## Index the genotype file 
index.genotype(genotype.f)
