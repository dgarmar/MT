##
## Plot sQTL sharing and tissue specificity results
##

## WD
setwd("/nfs/users/rg/dgarrido/sharing/blueprint/results")

## load pi1 values for all tissues
pi1.fdr1 <- as.matrix(read.table("pi1.1fdr.txt",sep="\t", header=TRUE))
pi1.fdr5 <- as.matrix(read.table("pi1.5fdr.txt",sep="\t", header=TRUE))

## Prepare for plotting
colnames(pi1.fdr1)=c("Monocytes","Neutrophils","T-cells (mRNA)","T-cells (total RNA)")
rownames(pi1.fdr1)=colnames(pi1.fdr1)
colnames(pi1.fdr5)=colnames(pi1.fdr1)
rownames(pi1.fdr5)=colnames(pi1.fdr5)
library(corrplot)

##Colors
col2 <- colorRampPalette(c("#67001F" ,"#D6604D","#D6604D","#F4A582", "#FDDBC7",
                           "#D1E5F0","#D1E5F0", "#053061"))  


## Plot sharing results
corrplot(t(pi1.fdr1),order="alphabet",method = "color",tl.col="black",cl.lim=c(0.7,1),
         col=col2(200),cl.length=6,is.corr=FALSE,addCoef.col=T)

corrplot(t(pi1.fdr5),order="alphabet",method = "color",tl.col="black",cl.lim=c(0.7,1),
         col=col2(200),cl.length=6,is.corr=FALSE,addCoef.col=T)
         
      

d<-function(m){
## Distance calculation as the correlation between the vectors of pi1 values corresponding to each two cell types
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

## Colors hash
hashCTCols<-new.env( );
hashCTCols[[ "Monocytes" ]]<-      "gold"
hashCTCols[[ "Neutrophils" ]]<- "darkseagreen"
hashCTCols[[ "T-cells (mRNA)"     ]]<-  "chocolate2"
hashCTCols[[ "T-cells (total RNA)"    ]]<-       "chocolate"

getPalette = colorRampPalette(c("grey20","grey96"))

colors=sapply(colnames(pi1.fdr5),function (x) hashTissueCols [[ as.character(x) ]]) #  (fdr 5)

## Plot heatmap and clustering
hv <- heatmap(t(pi1.fdr5),ColSideColors = colors, RowSideColors = colors,distfun=d,col=rev(getPalette(256)), margins = c(5,10),symm=T)


## Plot cell type specificity
ctspec.1fdr <- apply(t(pi1.fdr1), 1, function(x) 1-mean(x))
ctspec.5fdr <- apply(t(pi1.fdr5), 1, function(x) 1-mean(x))

df<-as.data.frame(ctspec.5fdr)


p <- ggplot(df, aes(rownames(df),df$ctspec.5fdr*100,fill=rownames(df)))+ geom_bar(stat="identity")
p + scale_fill_manual(breaks=colnames(df),values=colors)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Cell types") + ylab("Cell type specificity ( 1-mean(pi1) )") + ggtitle("Cell type specificity of sQTLs at 5%FDR")+
  theme(plot.title = element_text(size=20, face="bold"))+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.4))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))
