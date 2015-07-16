
##
## Heatmaps and clustering of tissues regarding their sharing pattern
##

d<-function(m){
 # Computes the distance between each two tissues as the correlation between the vectors of pi1 values corresponding
 # to each of them, estimates of the degree of sharing.
  h <- matrix(0,dim(m)[1],dim(m)[2])
  colnames(h)<-colnames(m)
  rownames(h)<-rownames(m)
  for (i in 1:dim(m)[1]){
    for (j in i:dim(m)[1]){
      h[j,i] <- 1-cor(m[i,],m[j,],method="pearson")
    }
  }
  return(as.dist(h))
}

## Hash of colors (color key established in GTEx)

hashTissueCols<-new.env( );
hashTissueCols[[ "Adipose" ]]<-      "tan1"
hashTissueCols[[ "Adrenal" ]]<- "darkseagreen"
hashTissueCols[[ "Whole"     ]]<-  "magenta"
hashTissueCols[[ "Artery"    ]]<-       "red"
hashTissueCols[[ "Brain"     ]]<-  "yellow2"
hashTissueCols[[ "Breast"   ]]<-      "cyan3"
hashTissueCols[[ "Colon"  ]]<-  "burlywood2"
hashTissueCols[[ "Cells" ]]<- "violet"
hashTissueCols[[ "Esophagus"  ]]<-  "burlywood3"
hashTissueCols[[ "Fallopian"  ]]<-  "mistyrose2"
hashTissueCols[[ "Fibroblasts" ]]<- "lightblue3"
hashTissueCols[[ "Heart" ]]<-"mediumorchid4"
hashTissueCols[[ "Kidney"   ]]<-    "bisque3"
hashTissueCols[[ "Liver"     ]]<-  "bisque3"
hashTissueCols[[ "Lung"   ]]<- "olivedrab3"
hashTissueCols[[ "Muscle"  ]]<-  "slateblue2"
hashTissueCols[[ "Minor"  ]]<-  "green"
hashTissueCols[[ "Nerve"   ]]<-      "gold1"
hashTissueCols[[ "Ovary"   ]]<-  "lightpink"
hashTissueCols[[ "Pancreas"  ]]<-  "goldenrod3"
hashTissueCols[[ "Pituitary" ]]<-"darkseagreen2"
hashTissueCols[[ "Prostate"  ]]<-      "grey85"
hashTissueCols[[ "Skin"   ]]<- "dodgerblue"
hashTissueCols[[ "Spleen"   ]]<- "bisque3"
hashTissueCols[[ "Small"   ]]<- "burlywood1"
hashTissueCols[[ "Stomach" ]]<-   "burlywood1"
hashTissueCols[[ "Testis"  ]]<-      "grey65"
hashTissueCols[[ "Thyroid" ]]<- "springgreen4"
hashTissueCols[[ "Uterus"  ]]<-  "mistyrose2"
hashTissueCols[[ "Vagina"  ]]<-  "mistyrose2"

## WD
setwd("/nfs/users/rg/dgarrido/sharing/results2")

## Get pi1 values for all tissues
pi1.fdr1<-as.matrix(read.table("pi1.1fdr.txt2",sep="\t", header=TRUE)) # txt3 -> 0.5, txt2 -> 1
rownames(pi1.fdr1)=colnames(pi1.fdr1)
pi1.fdr5<-as.matrix(read.table("pi1.5fdr.txt2",sep="\t", header=TRUE))
rownames(pi1.fdr5)=colnames(pi1.fdr5)

names=gsub("[ ._].*","",perl=TRUE,colnames(pi1.fdr5))
colors=sapply(names,function (x) hashTissueCols [[ as.character(x) ]])

require(graphics); require(grDevices)
colnames(pi1.fdr5)<-gsub("[-_.]"," ",perl=TRUE,colnames(pi1.fdr5))
rownames(pi1.fdr5)<-colnames(pi1.fdr5)

## Heatmap
hv <- heatmap(t(pi1.fdr5),ColSideColors = colors, RowSideColors = colors,distfun=d, col = cm.colors(256),margins = c(17,12),symm=T,cexRow=1.1,cexCol=1.1)

## For tissues with more than 100 samples
tissues.100 = sort(c("Muscle - Skeletal","Whole Blood","Skin - Sun Exposed (Lower leg)","Adipose - Subcutaneous","Artery - Tibial","Lung","Thyroid","Cells - Transformed fibroblasts","Nerve - Tibial","Esophagus - Mucosa","Esophagus - Muscularis","Artery - Aorta","Skin - Not Sun Exposed (Suprapubic)","Heart - Left Ventricle","Adipose - Visceral (Omentum)","Breast - Mammary Tissue","Stomach","Colon - Transverse","Heart - Atrial Appendage","Testis","Pancreas","Esophagus - Gastroesophageal Junction","Adrenal Gland","Colon - Sigmoid","Artery - Coronary","Cells - EBV-transformed lymphocytes","Brain - Cerebellum","Brain - Caudate (basal ganglia)","Liver"))
tissues.100 <- gsub(" - ",".",tissues.100)
tissues.100 <- gsub("-",".",tissues.100)
tissues.100 <- chartr(" ", "_",tissues.100)
tissues.100 <- gsub("[ _]\\(.+\\)","",tissues.100, perl=TRUE)

## Get pi1 values for tissues with more than 100 samples
pi1.1fdr.100 <- pi1.fdr1[tissues.100,tissues.100]
colnames(pi1.1fdr.100)<-gsub("[._]"," ",perl=TRUE,colnames(pi1.1fdr.100))
rownames(pi1.1fdr.100)<-gsub("[._]"," ",perl=TRUE,colnames(pi1.1fdr.100))

pi1.5fdr.100 <- pi1.fdr5[tissues.100,tissues.100]
colnames(pi1.5fdr.100)<-gsub("[._]"," ",perl=TRUE,colnames(pi1.5fdr.100))
rownames(pi1.5fdr.100)<-gsub("[._]"," ",perl=TRUE,colnames(pi1.5fdr.100))

names2=gsub("[ ._].*","",perl=TRUE,colnames(pi1.1fdr.100))
colors2=sapply(names2,function (x) hashTissueCols [[ as.character(x) ]])

## Heatmap
hv <- heatmap(t(pi1.5fdr.100),ColSideColors = colors2, RowSideColors = colors2,distfun=d, col = cm.colors(256), margins = c(5,10),symm=T)
