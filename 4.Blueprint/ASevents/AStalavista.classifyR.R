## Using AStalavista codes and directly from sQTLseekeR output, obtains the proportion of sQTLs associated to
## each event.

args = commandArgs(TRUE)

input.file = args[1] 
## Input tab-delimited file with columns 'transId1' and 'transId2', 
## the ids of the two transcript for which we want to compare the structure

tr.file = args[2]
## Input file with the exonic structure of the transcripts, created from 
## the annotation in previous steps ('getTransStructFromAnno.pl')

out.file = args[3] 
## The output file: 'input.file' with new columns informing about 
## the splicing events

out.stats = args[4] 
## Output file with a summary of which events were found.

trCoord = read.table(tr.file,header=TRUE,as.is=TRUE,fill=TRUE)
rownames(trCoord) = trCoord$txId

## Define functions to...

## 1.Derive splicing events code from two sets of exon (or UTR) boundary position
find.events <- function(tr1.pos,tr2.pos,check.5p=FALSE,check.3p=FALSE,posStrand=TRUE){
  pos.df = data.frame(pos=sort(unique(c(tr1.pos,tr2.pos)),decreasing=!posStrand))
  pos.df$tr1 = pos.df$pos %in% tr1.pos
  pos.df$tr2 = pos.df$pos %in% tr2.pos
  ii.3p = c(max(c(which(pos.df$tr1),1)),max(c(which(pos.df$tr2),1)))
  ii.5p = c(min(c(which(pos.df$tr1),1)),min(c(which(pos.df$tr2)),1))
  sep.c = c("-","^")
  s.1 = s.2 = 1
  c.cpt = 1
  c.tr1 = c.tr2 = ""
  ev.l = NULL
  for(ii in 1:nrow(pos.df)){
    if(pos.df[ii,"tr1"]!=pos.df[ii,"tr2"]){
      if(pos.df[ii,"tr1"]){
        c.tr1 = paste(c.tr1,c.cpt,sep.c[s.1],sep="")
      } else {
        c.tr2 = paste(c.tr2,c.cpt,sep.c[s.2],sep="")
      }
      c.cpt = c.cpt + 1
      if(check.5p & ii.5p[1]==ii) c.tr1 = paste("(",c.tr1,sep="")
      if(check.5p & ii.5p[2]==ii) c.tr2 = paste("(",c.tr2,sep="")
      if(check.3p & ii.3p[1]==ii) c.tr1 = paste(c.tr1,")",sep="")
      if(check.3p & ii.3p[2]==ii) c.tr2 = paste(c.tr2,")",sep="")
    } else if(c.cpt != 1){
      ev.l = c(ev.l, paste(c.tr1,c.tr2,sep=","))
      c.cpt = 1
      c.tr1 = c.tr2 = ""
    }
    if(pos.df[ii,"tr1"]) s.1 = 3 - s.1
    if(pos.df[ii,"tr2"]) s.2 = 3 - s.2
  }
  if(c.cpt != 1) ev.l = c(ev.l, paste(c.tr1,c.tr2,sep=","))
  ev.l
}

## 2.Translate a splicing code into a splicing event names
translate.event <- function(events){
  ev.tr = c(
    ",1-2^"="exon skipping",
    ",1^2-"="intron retention",
    "1-2^,3-4^"="mutually exclusive exon",
    "1-,2-"="alt 5'",
    "1^,2^"="alt 3'",
    "<>,(1-2^,(3-4^"="alt 5' UTR",
    "<>,1-2^),3-4^)"="alt 3' UTR",
    "<>,1^),2^)"="tandem 3' UTR",
    "<>,(1-,(2-"="tandem 5' UTR",
    "(1-2^,(3-4^"="alt first exon",
    "1-2^),3-4^)"="alt last exon",
    "1-2^,"="exon skipping",
    "1^2-,"="intron retention",
    "3-4^,1-2^"="mutually exclusive exon",
    "2-,1-"="alt 5'",
    "2^,1^"="alt 3'",
    "<>,(3-4^,(1-2^"="alt 5' UTR",
    "<>,3-4^)1-2^)"="alt 3' UTR",
    "<>,2^),1^)"="tandem 3' UTR",
    "<>,(2-,(1-"="tandem 5' UTR",
    "(3-4^(1-2^"="alt first exon",
    "3-4^),1-2^)"="alt last exon")
  ev.res = rep("complex event",length(events))
  ev.res[grepl("\\(",events)] = "complex event 5'"
  ev.res[grepl("\\)",events)] = "complex event 3'"
  ev.res[events %in% names(ev.tr)] = ev.tr[events[events %in% names(ev.tr)]]
  ev.res
}

## 3.Retrieve splicing code and event names from two transcript ids
classify.event <- function(tr1 = "ENST00000451283.1",tr2 = "ENST00000215882.5"){
  if(!any(tr1 %in% rownames(trCoord)) | !any(tr2 %in% rownames(trCoord)))
    return(list(class.df=data.frame(),class.stats=character(0)))
  tr1.cds.pos = as.integer(unlist(strsplit(c(trCoord[tr1,"cdsStarts"],trCoord[tr1,"cdsEnds"]),",")))
  tr2.cds.pos = as.integer(unlist(strsplit(c(trCoord[tr2,"cdsStarts"],trCoord[tr2,"cdsEnds"]),",")))
  mean.mm <- function(e)mean(c(min(e),max(e)))
  tr.center = c(mean.mm(tr1.cds.pos),mean.mm(tr2.cds.pos))
  posStrand = trCoord[tr1,"strand"]=="+"
  tr1.utr.pos = as.integer(unlist(strsplit(c(trCoord[tr1,"utrStarts"],trCoord[tr1,"utrEnds"]),",")))
  tr2.utr.pos = as.integer(unlist(strsplit(c(trCoord[tr2,"utrStarts"],trCoord[tr2,"utrEnds"]),",")))
  tr1.utr5.pos = tr1.utr.pos[ifelse(posStrand,1,-1)*(tr1.utr.pos-tr.center[1])<0]
  tr2.utr5.pos = tr2.utr.pos[ifelse(posStrand,1,-1)*(tr2.utr.pos-tr.center[2])<0]
  tr1.utr3.pos = tr1.utr.pos[ifelse(posStrand,1,-1)*(tr1.utr.pos-tr.center[1])>0]
  tr2.utr3.pos = tr2.utr.pos[ifelse(posStrand,1,-1)*(tr2.utr.pos-tr.center[2])>0]
  ev = find.events(tr1.cds.pos,tr2.cds.pos,check.5p=TRUE,check.3p=TRUE,posStrand=posStrand)
  if(length(c(tr1.utr5.pos,tr2.utr5.pos))>0)
    ev = c(ev,paste("<>",find.events(tr1.utr5.pos,tr2.utr5.pos,check.5p=TRUE,posStrand=posStrand),sep=","))
  if(length(c(tr1.utr3.pos,tr2.utr3.pos))>0)
    ev= c(ev, paste("<>",find.events(tr1.utr3.pos,tr2.utr3.pos,check.3p=TRUE,posStrand=posStrand),sep=","))
  ev.t = unique(translate.event(ev))
  list(class.df = data.frame(trPair=paste(tr1,tr2,sep="-"),classCode=paste(unique(ev),collapse=";"),classEvent=paste(sort(ev.t),collapse=";")),
       class.stats = ev.t, anySpl=!all(grepl("\\(",ev) | grepl("\\)",ev)))
}

## Main

res.in = read.table(input.file,header=TRUE,as.is=TRUE)
tr.to.class = unique(t(apply(res.in[,c("tr.first","tr.second")],1,sort)))
class.res = apply(tr.to.class,1,function(trs)classify.event(trs[1],trs[2]))
library(plyr)
class.df = ldply(class.res,function(l)l$class.df)
class.stats = ldply(class.res,function(l)if(length(l$class.stats)>0) data.frame(count=1,event=l$class.stats))
class.stats = aggregate(count~event,data=class.stats,sum)
class.stats$prop = class.stats$count / sum(class.stats$count)
class.stats$prop.sqtl = class.stats$count / nrow(tr.to.class)
spl.count = sum(unlist(lapply(class.res,function(l)if(length(l$class.stats)>0) l$anySpl)))
class.stats = rbind(class.stats, data.frame(event="any splicing",count=spl.count,prop=NA,prop.sqtl=spl.count/nrow(tr.to.class)))

res.in$trPair = apply(res.in[,c("tr.first","tr.second")],1,function(trs)paste(sort(trs),collapse="-"))
res.in = merge(res.in,class.df,all.x=TRUE)
res.in$trPair = NULL
write.table(res.in,file=out.file,quote=FALSE,row.names=FALSE,sep="\t")
write.table(class.stats,file=out.stats,quote=FALSE,row.names=FALSE,sep="\t")
