##
## Plot closeness
##

closeness <- function(file){
  # Cumulative proportion of intronic sQTLs according to their closeness to the closest exon
  
  table<-read.table(file)
  distances<-as.factor(table[,2])
  df<-as.data.frame(table(factor(distances)))
  df<-df[-c(1),]
  df<-rbind(c(0,0),df)
  
  total= sum(df$Freq)
  df$h = df$Freq/total
  
  for (i in 1:nrow(df)){
    if (i==1){
      Acum <-df[i,3]
      v<-c(Acum)
    }else{
      Acum <- Acum + df[i,3]
      v<-c(v,Acum)
    }
  }
  
  df<-cbind(df,v)
  return(df)
}

## Arguments from the command line
args <- commandArgs(TRUE)
sqtls.df<-closeness(args[1]) # closeness.txt
control.df<-closeness(args[2]) # closeness.control.txt

## PNG plot
png(args[3])
plot(type="l",as.character(sqtls.df[1:5001,1]),sqtls.df[1:5001,4],col="darkgreen",ylim=c(0, 1), xlim=c(0,5000), main="Intronic SNPs", xlab="Distance to the closest exon (bp)",ylab="Cumulative proportion")
lines(type="l",control.df[1:5001,4],col="darkred")
dev.off()
