
##
## Write a file with sqtlseeker.R results: all the associations gene-SNP and
## their p-values
##


 wd <- getwd()          # Save the WD from which the script is run
 setwd("../../../bin")  # Change WD. BatchJobs needs to be loaded in the path 
                        # where the configuration files are located

## Look at the paths in the specified order when loading a library  
.libPaths(c("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1",
            "/nfs/software/R/packages"))
library(BatchJobs)
library(sQTLseekeR)

setwd(wd) # Back to the initial WD

## Write a file with all the associations and p-values 
## (first check sQTL directory existence)

if (!file.exists("sQTL")){
 cat("## WARNING: sQTL directory is not available for this tissue. Any result can be reduced.")
 cat("\n")
 quit("no")
}

sQTL.reg <- makeRegistry(id="sQTL", seed=123, file.dir="sQTL")
res.f = "sQTLs-all.tsv"
if(file.exists(res.f)) file.remove(res.f)
tmp = reduceResultsList(sQTL.reg, fun=function(job, res){
  write.table(res, file=res.f, quote=FALSE, row.names=FALSE, col.names=!file.exists(res.f), append=file.exists(res.f), sep="\t")
})
