
## Needs to be present in the working directory when BatchJobs package is loaded. 
## It loads another R script file (here makeClusterFunctionsAdaptive.R) with the
## functions to use, and defines the cluster template.

 source("makeClusterFunctionsAdaptive.R")
 cluster.functions <- makeClusterFunctionsAdaptive("cluster.tmpl")
 mail.start <- "none"
 mail.done <- "none"
 mail.error <- "none"
 mail.from <- "<email@crg.eu>"
 mail.to <- "<email@crg.eu>"

