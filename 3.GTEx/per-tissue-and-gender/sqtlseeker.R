
##
##  Association test between genotypes and splicing ratios.
##  Called from multi_sqtlseeker.R, it uses sQTLseekeR's sqtl.seeker function to
##  perform the statistical test for association between each pair gene - SNP.
##  It uses R package BatchJobs to send the individual jobs to the computing 
##  cluster.
##

 wd <- getwd()      # Save the WD from which the script is run
 setwd("../../bin") # Change WD. BatchJobs needs to be loaded in the path 
                   # where the configuration files are located

## Look at the paths in the specified order when loading a library

 .libPaths(c("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1",
            "/nfs/software/R/packages"))
 library(BatchJobs)
 library(sQTLseekeR)

 setwd(wd)  # Back to the initial WD

## Input files: gene location and genotype information in required formats

 gene.bed.f <- "/nfs/users/rg/dgarrido/run-GTEx/input/genes.bed"
 genotype.indexed.f <- "/nfs/users/rg/dgarrido/run-GTEx/input/snps.tsv.bgz" 
 ## The genotype file is already indexed

## Test for gene/SNP associations

 # system("rm -rf sQTL") ## For debugging. Run to clean previous computations

## Make registry
 sQTL.reg <- makeRegistry(id="sQTL", seed=123, file.dir="sQTL")

## Load gene annotation file (BED)
 gene.bed <- read.table(gene.bed.f, as.is=TRUE, sep="\t")
 colnames(gene.bed) <- c("chr", "start", "end", "geneId")

## Load the result of prepare.trans.exp.R (splicing ratios of transcripts 
## suitable for the analysis) from the previously created registry 
 prepTE.reg <- makeRegistry(id="prepTE", seed=123, file.dir="prepTE")
 tre.df <- loadResult(prepTE.reg, 1)
 colnames(tre.df)[3:ncol(tre.df)] <- gsub(pattern="(.*)", # Adjust colnames
                                        replacement="GTEX\\1", perl=TRUE,
                                        colnames(tre.df)[3:ncol(tre.df)]) 

## Check
 if (dim(tre.df)[1]==0){
    cat("## WARNING: Skipped sample: sQTLseekeR cannot be run in this sample, 
      0 transcripts over the treshold in prepare.trans.exp")
    quit("no")
 }

## Subset genes (select the ones present in tre.df) and make chunks of 30 genes
 gene.bed <- subset(gene.bed, geneId %in% tre.df$geneId)
 nb.gene.per.chunk <- 30
 gene.chunks <- tapply(gene.bed$geneId, 
                     rep(1:ceiling(nrow(gene.bed)/nb.gene.per.chunk),
                         each=nb.gene.per.chunk)[1:nrow(gene.bed)], identity)

## Function that will be mapped
 sQTL.f <- function(chunk.id, imF){
    ## Subset transcripts for genes in a chunk and run sqtl.seeker function
    
    load(imF)
    genes <- gene.chunks[[chunk.id]]
    tre.df <- subset(tre.df, geneId %in% genes)
    gene.bed <- subset(gene.bed, geneId %in% genes)
    
    ## Look at the paths in the specified order when loading a library
    .libPaths(c("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1",
                "/nfs/software/R/packages"))  
    library(sQTLseekeR)
    
    ## Run
    sqtl.seeker(tre.df, genotype.indexed.f, gene.bed, genic.window=5e3, 
                svQTL=TRUE)
 }

 imF <- "sQTL-BJ-temp.RData" # Temporary file to store the datasets needed
 save(gene.chunks, tre.df, genotype.indexed.f, gene.bed, file=imF)

## Apply (map) the previous function to the created registry, and send the job 
## to the computing cluster, specifying the necessary resources.
 
 batchMap(sQTL.reg, sQTL.f, 1:length(gene.chunks), more.args=list(imF=imF))
 submitJobs(sQTL.reg, findNotDone(sQTL.reg), 
           resources=list(walltime="6:0:0", mem="4G", cores="1",
                          queue="short-sl65,long-sl65,rg-el6"), 
           wait=function(retries) 100, max.retries=10)


 ## Each chunk will constitute a job sent independently to the cluster.
 ## All the jobs will run in a parallel manner.
 
 # showStatus(sQTL.reg) 
 
