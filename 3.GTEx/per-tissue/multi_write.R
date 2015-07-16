
##
##  Iterate through tissue folders to run write.table.R in a per tissue manner
##

## WD
setwd("/nfs/users/rg/dgarrido/run-GTEx/")

## Getting tissue (folder) names 
sample.tissues <- read.table(file="input/sample-groups.txt", 
                             header=TRUE, sep="\t", as.is=TRUE)
tissues <- gsub(" - ", "-", sort(unique(sample.tissues[,3])[-55]))
tissues <- chartr(" ", "_", tissues)
tissues <- gsub("[ _]\\(.+\\)", "", tissues, perl=TRUE)

## Iterate and run write.table.R for each tissue
for (tissue in tissues){
  setwd(sprintf("results/%s", tissue))
  system("Rscript ../../bin/write.table.R 2>> stderr.txt >> stdout.txt")
  setwd("../..")
}


