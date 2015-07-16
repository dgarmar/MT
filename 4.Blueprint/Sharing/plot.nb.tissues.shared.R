

##
## Plot the number of cell types among which sQTLs are indeed shared
##

## WD
setwd("/nfs/users/rg/dgarrido/run-blueprint/run/results/")  # 5FDR, cufflinks

## Load shared sQTL counts
mon <- read.table("monocytes/cufflinks/sharing/sharing.counts.5fdr",as.is=T,header=F)
neut <- read.table("neutrophils/cufflinks/sharing/sharing.counts.5fdr",as.is=T,header=F)
tcm <- read.table("tcells-mrna/cufflinks/sharing/sharing.counts.5fdr",as.is=T,header=F)
tct <- read.table("tcells-totalrna/cufflinks/sharing/sharing.counts.5fdr",as.is=T,header=F)


#Prepare for plotting

names<- c("Monocytes","Neutrophils","T-Cells(mRNA)","T-Cells")
props <- data.frame(NULL)

j<-1
for (celltype in list(mon, neut, tcm, tct)){
  
  #celltype <- celltype[celltype[,1]!=1,]
  p <- data.frame(table(celltype[,1])/sum(table(celltype[,1])))
  #print(p)
  p <- cbind(rep(names[j],3),p)
  props <- rbind(props,p)
  j<-j+1
  
}

colnames(props)<-c("celltype","number","prop")

library(ggplot2)

## Plot

ggplot(props, aes(celltype, prop, fill=number)) + 
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.4))+
  scale_fill_brewer(palette = 17)+ 
  ggtitle("Pattern of effective sGene sharing")+
  theme(plot.title = element_text(size=20, face="bold"))+
  xlab("Percentage of sGenes")+
  ylab("Cell types")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.title = element_text(size=15, 
                                    face="bold"))+
  labs(fill = "Number of tissues")
