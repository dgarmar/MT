
##
##  Iterate through tissue folders to run write.table.R in a per tissue and per gender manner
##

## WD
setwd("/nfs/users/rg/dgarrido/run-GTEx/")

## Getting tissue (folder) names 
sample.tissues <- read.table(file="input/sample-groups.txt", 
                      header=TRUE,sep="\t", as.is=TRUE)
tissues <- gsub(" - ","-",sort(unique(sample.tissues[,3])))
tissues <- chartr(" ", "_",tissues)
tissues <- gsub("[ _]\\(.+\\)","",tissues, perl=TRUE)


## Iterate and run write.table.R for each tissue
for (tissue in tissues){
  setwd(sprintf("results/%s/female",tissue))
  system("Rscript ../../../bin/write.table.R 2>> stderr.txt >> stdout.txt")
  setwd("../male")
  system("Rscript ../../../bin/write.table.R 2>> stderr.txt >> stdout.txt")
  setwd("../../..")
}
