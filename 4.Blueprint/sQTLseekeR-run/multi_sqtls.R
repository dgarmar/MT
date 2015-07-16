
##
##  Iterate through cell type folders to run prepare.trans.exp.R in a 
##  per cell type and per quantifier(Flux,Cufflinks) manner 
##

## WD
setwd("/nfs/users/rg/dgarrido/run-blueprint/run/")

## Getting cell type (folder) names 
sample.tissues <- read.table(file="input/samples-celltypes.txt", header=TRUE,sep="\t", as.is=TRUE)
celltypes <- sort(unique(sample.tissues[,2]))

## Iterate and run sqtlseeker.R for each cell type and quantifier
for(celltype in celltypes){
  setwd(sprintf("results/%s/flux",celltype))
  system("qsub ../../../bin/sqtls.sh 2>> stderr.txt >> stdout.txt")
  setwd("../cufflinks")
  system("qsub ../../../bin/sqtls.sh 2>> stderr.txt >> stdout.txt")
  setwd("../../..")
}
