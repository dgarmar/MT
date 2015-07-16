##
##  Automatic plots for splicing QTLs
##  Can be used together with super.plotter.sh or manually selecting the cell type of interest
##

## args <- commandArgs(TRUE) # When super.plotter.sh is used

## Load BatchJobs to get the tre.df files that are in prepTE folders
setwd("/nfs/users/rg/dgarrido/run-blueprint/run/bin")
library(BatchJobs)
options(scipen=999) # Avoid scientific notation

## Load required files for plotting: 

load("../input/snps.samples.RData")
snps.safe = snps.samples

setwd("../results/tcells-totalrna/cufflinks") # The cell type you want to explore 
# OR setwd(sprintf("../results/%s",args[1])) when using super.plotter.sh

prepTE.reg = makeRegistry(id="prepTE", seed=123, file.dir="prepTE")
tre.df = loadResult(prepTE.reg, 1)

sig = read.table("sQTLs-sig5.tsv",header=TRUE,as.is=TRUE,sep="\t" ) #1%FDR

## Exploratory analysis
plot(sig$md,-log10(sig$qv),col="blue",pch=20) # Initial semi-volcano plot
sig.sorted<-sig[order(sig$md,-rank(sig$qv),decreasing=TRUE),] 
  #Sort by fold md value (decreasing) and p-val (increasing)

## Plot pdf

pdf("autoplot.pdf", paper = "a4r", width = 0, height = 0)

for (x in 1:1277){
  gene = sig.sorted[x,1]
  snp = sig.sorted[x,2]
  cat(sprintf("%s. %s - %s",x,gene,snp))
  cat("\n")
  
  subset.tr = tre.df[which(tre.df$geneId==gene),]
  subset.tr <- subset.tr[,colSums(is.na(subset.tr))<nrow(subset.tr)] # Remove columns with all NAs
  
  individuals = colnames(subset.tr)[-c(1:2)]
  subset.snp = as.numeric(snps.samples[which(snps.samples$snpId==snp),individuals])
  names(subset.snp) = individuals
  
  table(subset.snp) # TO CHECK
  
  rownames(subset.tr) = subset.tr$trId
  
  
  #   if (length(table(subset.snp))== 4){
  #     if (sum(table(subset.snp)[2:4]> sum(table(subset.snp)[2:4])/10) < 3){
  #       next
  #     }
  #   }
  
  AA = data.frame(t(subset.tr[,names(which(subset.snp==0))]))
  AB = data.frame(t(subset.tr[,names(which(subset.snp==1))]))
  BB = data.frame(t(subset.tr[,names(which(subset.snp==2))]))
  
  #   S = dim(AA)[1]+dim(AB)[1]+dim(BB)[1]
  #   
  #   if (length(table(subset.snp)) == 3){
  #     if (sum(table(subset.snp)) == S){
  #       if (sum(table(subset.snp) > S/10)< 3){
  #         next
  #       }
  #     }else{
  #       if (sum(table(subset.snp)[2:3] > S/10)< 2){
  #         next
  #       }
  #     }
  #   }else if (length(table(subset.snp)) == 2){
  #     if (sum(table(subset.snp) > S/10)< 2){
  #       next
  #     }
  #   }
  
  
  if (dim(AA)[1]==0){
    L = list(AB,BB)
    g = c("AB","BB")
    j = c(dim(AB)[1],dim(BB)[1])
    
  }else if (dim(AB)[1]==0){
    L=list(AA,BB)
    g=c("AA","BB")
    j = c(dim(AA)[1],dim(BB)[1])
    
  }else if (dim(BB)[1]==0){
    L=list(AA,AB)
    g=c("AA","AB")
    j = c(dim(AA)[1],dim(AB)[1])
    
  }else{
    L=list(AA,AB,BB)
    g=c("AA","AB","BB")
    j = c(dim(AA)[1],dim(AB)[1],dim(BB)[1])
  }
  
  
  library(ggplot2)
  library(gridExtra)
  
  i=1
  p <- list()
  for (geno in L){
    
    #   if(i==1) {g = "AA"}
    #   else if (i==2) {g = "AB"}
    #   else {g = "BB"}
    
    
    df<-data.frame()
    
    for (cn in colnames(geno)) {
      
      #     if( mean ( geno[,cn], na.rm=TRUE ) < 0.01 ){
      #       df.tmp <-data.frame(tr=rep("Others",nrow(geno)),y=geno[,cn])
      
      #     }else{
      df.tmp <-data.frame(tr=rep(cn,nrow(geno)),y=geno[,cn])
      #     }
      df<-rbind(df,df.tmp)
    }
    
    p[[i]] <- ggplot(df, aes(tr,y))+
      geom_boxplot(aes(fill=tr)) + labs(x = "Isoform") + labs(y = "relative abundance") +
      theme(legend.position="none")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      scale_x_discrete(labels=levels(df$tr))+ ggtitle(sprintf("%s(%s)",g[i],j[i]))+
      stat_summary(fun.y = "mean", geom = "point", pch= 20, size= 3, color= "white")+
      ylim(0,1)
    
    i=i+1
  }
  
  if (i==3){
    suppressWarnings(grid.arrange(p[[1]],p[[2]],ncol=2,main=sprintf("%s - %s",gene,snp)))
  }else if (i==4){
    suppressWarnings(grid.arrange(p[[1]],p[[2]],p[[3]],ncol=3,main=sprintf("%s - %s",gene,snp)))
    
  }
  
}

dev.off()
