
##
##  Iterate through cell type folders to run prepare.trans.exp.R in a 
##  per cell type and per quantifier(Flux,Cufflinks) manner 
##

## WD
setwd("/nfs/users/rg/dgarrido/run-blueprint/run/")

## Getting cell types
sample.celltypes <- read.table(file="input/samples-celltypes.txt", header=TRUE,sep="\t", as.is=TRUE)
celltypes <- sort(unique(sample.celltypes[,2]))

## Iterate and run prepare.trans.exp.R for each cell type and quantifier

for (celltype in celltypes){
  setwd(sprintf("results/%s/flux",celltype))
  system("Rscript ../../../bin/write.table.R 2>> stderr.txt >> stdout.txt")
  setwd("../cufflinks")
  system("Rscript ../../../bin/write.table.R 2>> stderr.txt >> stdout.txt")
  setwd("../../..")
}
