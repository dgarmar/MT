
##
##  Transcript quantification data pre-processing and obtention of splicing 
##  ratios. Called from multi_prepare.trans.exp.R. It uses sQTLseekeR's 
##  prepare.trans.exp function to filter out the genes and transcripts not 
##  suitable for the analysis. It uses BatchJobs R package to send individual
##  jobs to the computing cluster. 
##

 wd <- getwd()       # Save the WD from which the script is run
 setwd("../../bin")  # Change WD. BatchJobs needs to be loaded in the path 
                    # where the configuration files are located

## Look at the paths in the specified order when loading a library 

 .libPaths(c("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1",
            "/nfs/software/R/packages"))
 library(BatchJobs)
 library(sQTLseekeR)

 setwd(wd) # Back to the initial WD

 args <- commandArgs(TRUE) # Command line input (from multi_prepare.trans.exp.R)
 t.index <- as.numeric(args[1])

## Input file
 trans.exp.f <- "../../input/te.df.RData" # Transcript quantification data file 
                                        # in the required format

## Getting the IDs of samples of the tissue of interest
## Note: This is relevant here because we will study different 
##       subsets of samples (corresponding to different tissues).

 sample.groups <- read.table(file="../../input/sample-groups.txt", 
                           header=TRUE, as.is = TRUE, sep="\t")
 tissue.sel <- sort(unique(sample.groups[,3])[-55])[t.index]
 subset.samples <- subset(sample.groups, tissue==tissue.sel[1])$sample

## Getting the names of samples in the genotype file

 snps.header <- read.table(file="../../input/snps.tsv.header.X",
                         as.is=TRUE, sep="\t")
 snps.header <- snps.header[-c(1:2),1] 
 snps.header <- gsub(pattern="GTEX-([^-]+)-.*", replacement="\\1", 
                   perl=TRUE, snps.header)

## Transcript quantification data pre-processing

 # system("rm -rf prepTE") ## For debugging. Run to clean previous computations

## Make registry
 prepTE.reg <- makeRegistry(id="prepTE", seed=123, file.dir="prepTE") 

## Define function that will be mapped
 prepTE.f <- function(te.file, samples, header){ 
    ## Subset samples and run prepare.trans.exp
    
    load(te.file)
    colnames(te.df)[1:2] <- c("trId", "geneId")
    
    ## Take the samples of the tissue of interest
    te.df <- te.df[,c(1,2,which(colnames(te.df)%in%samples))] 
    
    colnames(te.df) <- gsub(pattern="GTEX-([^-]+)-.*", replacement="\\1", 
                            perl=TRUE, colnames(te.df)) # Adjust colnames
    
    ## Take the samples with associated genotype
    te.df <- te.df[,c(1,2,which(colnames(te.df)%in%header))] 
    
    ## Check
    if(!identical(sum(!is.na(colnames(te.df))),
                  sum(!is.na(unique(colnames(te.df)))))){
        cat("Error occurred: repeated colnames")
        quit("no")
    }   
    
    ## Run
    .libPaths(c("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1",
                "/nfs/software/R/packages"))
    library(sQTLseekeR)
    prepare.trans.exp(te.df, min.transcript.exp=0.01, min.gene.exp=0.01, 
                      min.dispersion=0.1)
}

## Apply (map) the previous function to the created registry, and send the job 
## to the computing cluster, specifying the necessary resources.

 batchMap(prepTE.reg, prepTE.f, trans.exp.f, 
         more.args=list(samples=subset.samples, header=snps.header))
 submitJobs(prepTE.reg, 1, resources=list(walltime="6:0:0", cores="1", mem="50G",
                                         queue="long-sl65,short-sl65,rg-el6"), 
           wait=function(retries) 100, max.retries=10)

## Take into account that prepare.trans.exp.R is run from 
## multi_prepare.trans.exp.R as many times as tissues are analysed . A job 
## corresponding to each tissue will be generated. Once the execution of 
## multi_prepare.trans.exp.R has finished, all the jobs will be running in a 
## parallel manner

# showStatus(sQTL.reg)
