
##
## Print number of tested protein coding genes, per tissue
##

## WD
setwd("/nfs/users/rg/dgarrido/run-GTEx1/bin")
library(BatchJobs)

## Tissue (folder) names retrieval
sample.tissues <- read.table(file="../input/sample-groups.txt", header=TRUE,sep="\t", as.is=TRUE)
tissues <- gsub(" - ","-",sort(unique(sample.tissues[,3])[-55]))
tissues <- chartr(" ", "_",tissues)
tissues <- gsub("[ _]\\(.+\\)","",tissues, perl=TRUE)

## Protein coding genes
pcg <- read.table(file="../../run-GTEx/input/genes.protcod.bed", header=TRUE,sep="\t", as.is=TRUE)

## Print
for (tissue in tissues){
  setwd(sprintf("../results/%s",tissue))
  prepTE.reg = makeRegistry(id="prepTE", seed=123, file.dir="prepTE")
  print(tissue)
  tre.df<-loadResult(prepTE.reg,1)
  print(length(which(unique(tre.df[,2])%in%pcg[,4])))
  setwd("../../bin")
}
