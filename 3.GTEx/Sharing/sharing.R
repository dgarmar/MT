#!/usr/bin/R

#Compute sharing (pi1 estimates) thanks to qvalue package

## Input from the command line (sharing.R called from sharing.sh)
args <- commandArgs(TRUE)
pvalues<-as.numeric(args)

## load qvalue package
.libPaths("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1") 
library(qvalue)

## pi1 estimation
pi1 = 1-qvalue(pvalues)$pi0
cat(pi1)
