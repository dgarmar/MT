
##
## Plot sex-specificity estimates
##

## WD
setwd("/nfs/users/rg/dgarrido/sharing/fe-male")

## Input pi1 values
ts<-read.table("tiss.spec",as.is=T,sep="\t")
colnames(ts)<-c("tissue","female","male","sort")
ts<-ts[order(ts[,4],decreasing=TRUE),][,1:3]

library(plyr)
library(reshape2)
library(ggplot2)

ts$male<--ts$male
ts$tissue<-factor(ts$tissue,levels=ts$tissue,labels=ts$tissue)

ts.ok<-melt(ts,value.name='tisspec',variable.name = 'Sex',id.vars="tissue")

## Plot pdf

p<-ggplot(ts.ok,aes(x=tissue,y=tisspec,fill=Sex))
p +  geom_bar(subset = .(Sex == "male"), stat ="identity")+
  geom_bar(subset = .(Sex == "female"), stat ="identity")+
  theme_bw()+
  ggtitle("Sex specificity of sQTLs")+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Sex specificity (1-p1)")+
  xlab("Tissues")+
  theme(axis.text=element_text(size=12),
               axis.title=element_text(size=16,face="bold"))+
  guides(fill=guide_legend(title="Sex"))+
  theme(legend.title = element_text(colour="black", size=14, face="bold"))+
  theme(legend.text = element_text(colour="black", size = 12))

