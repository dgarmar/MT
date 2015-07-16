
##
##  Iterate through tissue folders to run sqtls.sh in a per tissue manner
##

## WD
setwd("/nfs/users/rg/dgarrido/run-GTEx/")

## Getting tissue (folder) names 
sample.tissues <- read.table(file="input/sample-groups.txt", 
                             header=TRUE,sep="\t", as.is=TRUE)
tissues <- gsub(" - ", "-", sort(unique(sample.tissues[,3])[-55]))
tissues <- chartr(" ", "_", tissues)
tissues <- gsub("[ _]\\(.+\\)", "", tissues, perl=TRUE)

## Iterate and run sqtls.sh for each tissue
for (tissue in tissues){
  setwd(sprintf("results/%s",tissue))
  system("qsub ../../bin/sqtls.sh 2>> stderr.txt >> stdout.txt")
  setwd("../..")
}
