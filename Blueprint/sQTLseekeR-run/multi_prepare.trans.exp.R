
##
##  Iterate through cell type folders to run prepare.trans.exp.R in a 
##  per cell type and per quantifier(Flux,Cufflinks) manner 
##

## WD
setwd("/nfs/users/rg/dgarrido/run-blueprint/run/")

## Getting cell types

sample.tissues <- read.table(file="input/samples-celltypes.txt",
                             header=TRUE,sep="\t", as.is=TRUE)
celltypes <- sort(unique(sample.tissues[,2]))


## Iterate and run prepare.trans.exp.R for each cell type

t.index=1
for (celltype in celltypes){
  system(sprintf("mkdir results/%s",celltype))
  setwd(sprintf("results/%s",celltype))
  for (q in 1:2){
    if (q==1){
      system("mkdir flux")
      setwd("flux")
    }else if (q==2){
      system("mkdir cufflinks")
      setwd("cufflinks")
    } 
    system(sprintf("Rscript ../../../bin/prepare.trans.exp.R %s %s 
                   2> stderr.txt > stdout.txt", t.index, q))
    setwd("..")
  }
  setwd("../..")
  t.index=t.index+1
}
