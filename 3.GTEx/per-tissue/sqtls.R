
##
## Perform FDR correction and get significant sQTLs. Called from multi_sqtls.R
##

wd <- getwd()       # Save the WD from which the script is run
setwd("../../bin")  # Change WD. BatchJobs needs to be loaded in the path 
                    # where the configuration files are located

## Look at the paths in the specified order when loading a library.
.libPaths(c("/nfs/users/rg/dgarrido/R/x86_64-redhat-linux-gnu-library/3.1",
            "/nfs/software/R/packages"))

library(BatchJobs)
library(sQTLseekeR)

setwd(wd) # Back to the initial WD

## Input file: result (all sQTLs)
res.df = read.table("sQTLs-all.tsv", header=TRUE, as.is=TRUE, sep="\t")

## Perform FDR correction and get significant sQTLs (in this case at 1% FDR)
 sig.f = "sQTLs-sig1.tsv"
 sqtls.df = suppressWarnings(sqtls(res.df, FDR=.01))
 write.table(sqtls.df, file=sig.f,quote=FALSE, row.names=FALSE,
            col.names=TRUE, sep="\t")
