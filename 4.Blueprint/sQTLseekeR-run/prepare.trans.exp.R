
##
##  Transcript quantification data pre-processing and obtention of splicing 
##  ratios. Called from multi_prepare.trans.exp.R. It uses sQTLseekeR's 
##  prepare.trans.exp function to filter out the genes and transcripts not 
##  suitable for the analysis. It uses BatchJobs R package to send individual
##  jobs to the computing cluster. 
##

 wd <- getwd()       # Save the WD from which the script is run
 setwd("../../..bin")  # Change WD. BatchJobs needs to be loaded in the path 
                     # where the configuration files are located

## Look at the paths in the specified order when loading a library 
.libPaths(c("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1",
                    "/nfs/software/R/packages"))

library(BatchJobs)
library(sQTLseekeR)


args <- commandArgs(TRUE) # Command line input (from multi_prepare.trans.exp.R)
setwd(wd) # Back to the initial WD
t.index <- as.numeric(args[1]) # Command line input (from multi_prepare.trans.exp.R:tissue)
q <- as.numeric(args[2])  # Command line input (from multi_prepare.trans.exp.R: gender)

# Input
if (q==1) {
  trans.exp.f = "~/run-blueprint/run/input/te.df.flux.RData"
}else if (q==2){
  trans.exp.f = "~/run-blueprint/run/input/te.df.cufflinks.RData"
}

## Getting the IDs of samples of celltype
## Note: This is relevant here because we will study different subset of samples (celltypes).
sample.groups = read.table(file="~/run-blueprint/run/input/samples-celltypes.txt", header=TRUE,
                    as.is = TRUE, sep="\t")
celltype.sel = sort(unique(sample.groups[,2]))[t.index]
subset.samples = subset(sample.groups,cell.type==celltype.sel)$sample

# Samples with genotype information
snps.header = as.character(read.table(file="~/run-blueprint/run/input/snps.tsv.header",
                                      as.is=TRUE,sep="\t"))[-c(1:4)]


## Transcript quantification data pre-processing

# system("rm -rf prepTE") ## For debugging. Run to clean previous computations

## Make registry
prepTE.reg <- makeRegistry(id="prepTE", seed=123, file.dir="prepTE")

## Define function that will be mapped
prepTE.f <- function(te.file, samples, header){
   ## Subset samples and run prepare.trans.exp
  load(te.file)
  
  # colnames(te.df)[1:2] = c("trId", "geneId") # Previously done
  ## Take the samples of the tissue of interest
  te.df = te.df[,c(1,2,which(colnames(te.df)%in%samples))] 
  
  ## Take the samples of with genotype
  colnames(te.df) <- gsub(pattern="(S.{5}).*", replacement="\\1", perl=TRUE,colnames(te.df)) 
  te.df <- te.df[,c(1,2,which(colnames(te.df)%in%header))] 
  
  ## Check
  if(!identical(length(colnames(te.df)),length(unique(colnames(te.df))))){
    cat("Error occurred: repeated colnames","\n")
    quit("no")
  }   
  
  ## Run
  .libPaths(c("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1","/nfs/software/R/packages"))
  library(sQTLseekeR)
  prepare.trans.exp(te.df)
}

## Apply (map) the previous function to the created registry, and send the job 
## to the computing cluster, specifying the necessary resources.
batchMap(prepTE.reg, prepTE.f,trans.exp.f, more.args=list(samples=subset.samples,header=snps.header))
submitJobs(prepTE.reg, 1, 
          resources=list(walltime="6:0:0", cores="1",mem="4G",queue="rg-el6,short-sl65,long-sl65"), 
          wait=function(retries) 100, max.retries=10)
          
## Take into account that prepare.trans.exp.R is run from 
## multi_prepare.trans.exp.R as many times as celltypes are analysed . A job 
## corresponding to each tissue will be generated. Once the execution of 
## multi_prepare.trans.exp.R has finished, all the jobs will be running in a 
## parallel manner

# showStatus(sQTL.reg)
