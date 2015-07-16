
##
##  Iterate through tissue folders to run prepare.trans.exp.R in a per tissue 
##  manner
##

## WD
 setwd("/nfs/users/rg/dgarrido/run-GTEx/")

## Getting tissue (folder) names 
 sample.tissues <- read.table (file="input/sample-groups.txt", 
                             header=TRUE, sep="\t", as.is=TRUE)
 tissues <- gsub(" - ", "-", sort(unique(sample.tissues[,3])[-55]))
 tissues <- chartr(" ", "_", tissues)
 tissues <- gsub("[ _]\\(.+\\)","", tissues, perl=TRUE)
 tissues.original <- sort(unique(sample.tissues[,3])[-55])

# system("rm -rf results/*") # For debugging

## Iterate and run prepare.trans.exp.R for each tissue
 
 t.index <- 1

 for (tissue in tissues){
    system(sprintf("mkdir results/%s", tissue))
    setwd(sprintf("results/%s", tissue))
    system(sprintf( "Rscript ../../bin/prepare.trans.exp.R %s 2> 
                 stderr.txt > stdout.txt", t.index ))
    setwd( "../.." )
    t.index <- t.index+1
 }

