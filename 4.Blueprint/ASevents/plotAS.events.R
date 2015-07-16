
# Plot proportions of sQTLs associated to each AS event

setwd("/nfs/users/rg/dgarrido/events/old/")

events <- read.delim("summary.1fdr.pc.nofilter")
events <- events[-dim(events)[1],c(1,3,4)]
events<-cbind(events,group=c(rep("sQTLseekeR",dim(events)[1])))

# Event names
names <-c("Tandem 5' UTR","Tandem 3' UTR","Mutually exclusive exons","Intron retention","Exon skipping",
          "Complex 5' event","Complex 3' event", "Complex event","Alternative last exon","Alternative first exon",
          "Alternative 5' UTR", "Alternative acceptor","Alternative 3' UTR", "Alternative donor")

# Colors
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9,"Blues"))(20)[7:20]
getPalette<-c("#75B3D8" ,"#62A8D2", "#4090C5", "#3282BD","#8BBFDC", "#519CCB",
              "#1966AD","#083D7F","#08306B", "#2474B6",  "#B0D2E7", "#A0CAE1","#0E59A2" , "#084B94" ) 

# Plot
library(RColorBrewer)
library(ggplot2)

setwd("~")
pdf("events.pdf", height=10, width=10)

ggplot(events, aes(event, prop.sqtl,fill=event)) + 
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Classification of major transcript changes")+
  theme(plot.title = element_text(size=19, face="bold"))+
  ylab("Proportion of sQTLs")+
  xlab("")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.title = element_text(size=15, 
                                    face="bold"))+
  coord_flip()+
  theme(axis.text.x = element_text(size=15,angle = 0, hjust = 0.4))+
  theme(axis.text.y = element_text(size=15))+
  scale_x_discrete(limits=rev(sort(unique(events$event))),labels=names)+
  scale_fill_manual(values = getPalette)+
  theme(legend.position="none")

dev.off()
