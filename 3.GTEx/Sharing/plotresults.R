
##
## Plot sQTL sharing and tissue specificity
##

## WD
setwd("/nfs/users/rg/dgarrido/sharing/results2")

## pi1 values for all tissues
pi1.fdr1<-as.matrix(read.table("pi1.1fdr.txt2",sep="\t", header=TRUE)) # txt3 -> 0.5, txt2 -> 1
pi1.fdr5<-as.matrix(read.table("pi1.5fdr.txt2",sep="\t", header=TRUE))

rownames(pi1.fdr1)=colnames(pi1.fdr1)
rownames(pi1.fdr5)=colnames(pi1.fdr5)

## pi1 values for tissues with more than 100 samples
tissues.100 = sort(c("Muscle - Skeletal","Whole Blood","Skin - Sun Exposed (Lower leg)","Adipose - Subcutaneous","Artery - Tibial","Lung","Thyroid","Cells - Transformed fibroblasts","Nerve - Tibial","Esophagus - Mucosa","Esophagus - Muscularis","Artery - Aorta","Skin - Not Sun Exposed (Suprapubic)","Heart - Left Ventricle","Adipose - Visceral (Omentum)","Breast - Mammary Tissue","Stomach","Colon - Transverse","Heart - Atrial Appendage","Testis","Pancreas","Esophagus - Gastroesophageal Junction","Adrenal Gland","Colon - Sigmoid","Artery - Coronary","Cells - EBV-transformed lymphocytes","Brain - Cerebellum","Brain - Caudate (basal ganglia)","Liver"))
tissues.100 <- gsub(" - ",".",tissues.100)
tissues.100 <- gsub("-",".",tissues.100)
tissues.100 <- chartr(" ", "_",tissues.100)
tissues.100 <- gsub("[ _]\\(.+\\)","",tissues.100, perl=TRUE)

pi1.1fdr.100 <- pi1.fdr1[tissues.100,tissues.100]
pi1.5fdr.100 <- pi1.fdr5[tissues.100,tissues.100]

colnames(pi1.1fdr.100)<-gsub("[._]"," ",perl=TRUE,colnames(pi1.1fdr.100))
colnames(pi1.5fdr.100)<-gsub("[._]"," ",perl=TRUE,colnames(pi1.1fdr.100))
rownames(pi1.1fdr.100)<-gsub("[._]"," ",perl=TRUE,colnames(pi1.1fdr.100))
rownames(pi1.5fdr.100)<-gsub("[._]"," ",perl=TRUE,colnames(pi1.1fdr.100))

## Color palettes
col1 <- colorRampPalette(c("#67001F","#67001F", "#B2182B","#B2182B", "#D6604D","#D6604D", "#F4A582", "#F4A582", "#FDDBC7", "#D1E5F0","#D1E5F0", "#053061"))	
col2 <- colorRampPalette(c("#67001F", "#B2182B","#B2182B", "#D6604D","#D6604D","#D6604D", "#F4A582","#F4A582", "#FDDBC7","#D1E5F0", "#053061"))  
col3 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7","#D1E5F0", "#053061"))  corrplot(pi1.1fdr.100, order="hclust",hclust.method="complete",method = "color", cl.lim=c(0.5,1),col=col2(200),tl.col="black",cl.length=10,is.corr=FALSE,mar = c(4,2,1,1))

## Heatmap
library(corrplot)
corrplot(t(pi1.1fdr.100),order="alphabet",method = "color", cl.lim=c(0.5,1),col=col3(200),tl.col="black",cl.length=10,is.corr=FALSE,mar = c(4,2,1,1))

## new WD
setwd("/nfs/users/rg/dgarrido/run-GTEx1/info")

# Load tissue abbreviated names
abbrev<-read.table("abbrev.txt",sep="\t", header=TRUE)
abbrev2<-abbrev[-c(7,23,25,26,32),] # Remove the 0 ones, bladder and so on
abbrev3<-abbrev[-c(7,8,9,11,13,14,15,16,17,18,19,20,23,25,26,32,35,38,41,43,44,47,48,52,53),]# Remove the <100 samples ones
nb.samples<-read.table("nb.samples",sep="\t", header=TRUE)

## Back to the old WD
setwd("/nfs/users/rg/dgarrido/sharing/results2")

## Compute tissue specificity
tisspec_all.1fdr <- apply(pi1.fdr1, 2, function(x) 1-mean(x))
tisspec_all.5fdr <- apply(pi1.fdr5, 2, function(x) 1-mean(x))  

## Plot (1%FDR) - all tissues
df<-as.data.frame(tisspec_all.1fdr)
df<-cbind(df,type=gsub("[ ._].*","",perl=TRUE,rownames(df)))
colors=sapply(unique(df$type),function (x) hashTissueCols [[ as.character(x) ]])
rownames(df)<-gsub("[._]"," ",perl=TRUE,rownames(df))

p <- ggplot(df, aes(rownames(df),df$tisspec_all.1fdr,colour=df$type)) + geom_point(cex=3)
p + scale_colour_manual(breaks=df$type,values=colors)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Tissues") + ylab("Tissue specificity ( 1-mean(pi1) )") + ggtitle("GTEx scaling-up - sQTLs Tissue specificity 1%FDR")

# Plot (5%FDR) - all tissues
df2<-as.data.frame(tisspec_all.5fdr)
df2<-cbind(df2,type=gsub("[ ._].*","",perl=TRUE,rownames(df2)))
colors2=sapply(unique(df2$type),function (x) hashTissueCols [[ as.character(x) ]])
rownames(df2)<-gsub("[._]"," ",perl=TRUE,rownames(df2))

p <- ggplot(df2, aes(rownames(df2),df2$tisspec_all.5fdr,colour=df2$type)) + geom_point(cex=3)
p + scale_colour_manual(breaks=df2$type,values=colors)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Tissues") + ylab("Tissue specificity ( 1-mean(pi1) )") + ggtitle("GTEx scaling-up - sQTLs Tissue specificity 5%FDR")

# Only >100 samples

tisspec.1fdr <- apply(pi1.1fdr.100, 2, function(x) 1-mean(x))
tisspec.5fdr <- apply(pi1.5fdr.100, 2, function(x) 1-mean(x))  

# Plot (1%FDR)
df3<-as.data.frame(tisspec.1fdr)
df3<-cbind(df3,type=gsub("[ ._].*","",perl=TRUE,rownames(df3)))
colors3=sapply(unique(df3$type),function (x) hashTissueCols [[ as.character(x) ]])

p <- ggplot(df3, aes(rownames(df3),df3$tisspec.1fdr,fill=df3$type)) + geom_bar(stat="identity")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme (axis.text.y = element_text(angle = 90, hjust = 0.5)) 
p + scale_fill_manual(breaks=df3$type,values=colors3)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Tissues") + ylab("Tissue specificity ( 1-mean(pi1) )") + ggtitle("GTEx scaling-up - sQTLs Tissue specificity 5%FDR")



