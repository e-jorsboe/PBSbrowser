
shinyDir<-commandArgs(trailingOnly=T)[1]
name<-commandArgs(trailingOnly=T)[2]

setwd(".")

theFile<-paste0(shinyDir,"/tmp.Rdata")

source("scripts/tmpFunctionsPBSV2.R") 

load(theFile)
attach(res)


## does the whole Genome PBS - takes a very long time, therefore put into seperate scripts
wgTable<-wholeGenomePBS(windows=windows,fall=fall,nInd=nInd,pos=pos,rs=rs,chr=chr,n=n,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs,pop1=pop1,pop2=pop2,pop3=pop3,winSize=winSize,minWin=minWin,shinyPBS=shinyPBS,SNPsinChr=SNPsinChr,maxChr=max(chr),FstOnly=(FstOnly=="YES"))
    
    
    
inInterval<-PBSTable(wgTable=wgTable,chr=chr,pop1=pop1,pop2=pop2,pop3=pop3,ifWindows=windows,winSize=winSize,minWin=minWin,genes=genes,all="YES",start=start,end=end,thisChr=11,maxChr=max(chr),FstOnly=(FstOnly=="YES"))
inInterval$pos<-as.numeric(inInterval$pos)
## might have to change with Fst
## inInterval[,c("PBS",paste("Fst",pop1,pop2,sep=""),paste("Fst",pop1,pop3,sep=""),paste("Fst",pop2,pop3,sep=""),"quantile")]<-round(inInterval[,c("PBS",paste("Fst",pop1,pop2,sep=""),paste("Fst",pop1,pop3,sep=""),paste("Fst",pop2,pop3,sep=""),"quantile")],5)


setwd(shinyDir)

con <- gzfile(sub("Rdata","gz",theFile),"w")
#con <- gzfile(name,"w")
write.table(inInterval,con,sep="\t",quote=F,row=F)
close(con)

system(paste("cp ",shinyDir,"/tmp.gz ",shinyDir,"/",name,sep=""))
