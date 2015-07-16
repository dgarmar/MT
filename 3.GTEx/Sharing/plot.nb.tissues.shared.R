
##
## Plot the number of tissues among which sQTLs are indeed shared
##

## WD
setwd("/nfs/users/rg/dgarrido/run-GTEx1/results.over100/")

## Load sQTL counts in tissues
sqtls.sig1 <- read.table("sqtls.sig1.vs.nb.tissues.over100.pc",as.is=T,header=F)
sqtls.sig5 <- read.table("sqtls.sig5.vs.nb.tissues.over100.pc",as.is=T,header=F)

## Prepare dataset
sqtls.sig1 <- sqtls.sig1[sqtls.sig1[,1]!=1,]
props <- data.frame(table(sqtls.sig1[,1])/sum(table(sqtls.sig1[,1])))
colnames(props)<-c("tissues","prop")

## Plot
library(ggplot2)
ggplot(props, aes(x=tissues,y=prop,fill=prop)) + 
  geom_bar(stat = "identity") + 
  scale_fill_gradient2(low='white', mid='lightskyblue1', high='dodgerblue3', space='Lab')+ 
  ggtitle("GTEx scale-up phase - sQTLs sharing")+
  theme(plot.title = element_text(size=20, face="bold"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.4))+
  ylab("Proportion of shared sQTLs")+
  xlab("Number of tissues")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  theme(legend.title=element_blank())
