passWord<-c("we<3Animals", "fiskfisk" )


require(rCharts)
options(RCHART_WIDTH = 800)

library(shiny)
source("/home/albrecht/Rfun/shiny.R")

setwd(".")

source("scripts/tmpFunctionsPBSV2.R") 

folderName<-scan("folderName",what="ADA")[1]

shinyDir<-"/pontus/data/shiny/emil/"
shinyPBS<-paste0(shinyDir,"pbs")
shinyCSV<-paste0(shinyDir,"pbs.csv")

read1 <- function(pop,popNames1,theFile,nSNP){
    
    w <- which(popNames1==pop)
    con<-file(theFile,"rb")
    seek(con,where=nSNP*(w-1)*8)
    r<-readBin(con,"double",nSNP,size=8)
    close(con)
    r
}

######################################
## read in data, freqs and individuals from admixture output
##############################

## nInd,nPop,nSNP
load("data/info.Rdata")
## chr,pos        
load("data/snpInfo.Rdata")

l<-list()
for(popName in popNames){
    
    l[[popName]]<-read1(popName,popNames,"data/freq.bin",nSNP)
    
}

fall<-do.call(cbind,l)
colnames(fall)<-popNames
names(nInd)<-popNames
SNPsinChr<-read.table("data/SNPinChr.txt",h=F,as.is=T,colClasses = "integer") ## number of SNPs in each chr used for whole genome win

## new genes file with only 4 columns chr, cdsStart, cdsEnd & name2
dat<-read.table("data/annoSelectedCols.txt",as.is=T,head=T)
geneFile<-"data/genes.anno"

########################################################################################

shinyServer(function(input, output, session) {
    
    
        PBSsinglePlot <- eventReactive(input$runPlot, {
              ##input <- list(start=60,end=62,chr=1,ifWindows="NO",winSize=50000,pop1="Ghana",pop2="Kenya",pop3="Uganda",minWin=10,together="NO",quantileLine=95,FstOnly="YES")
            
            if(pw(input,passWord)){
              return()
            }
            if(input$pop1==input$pop2 | input$pop1==input$pop3 | input$pop2==input$pop3 ){
              return(plot(1,main="You have picked on population at least twice!"))
            }
            start<-input$start*1e06
            end<-input$end*1e06
            
            ## clearing saved shiny saved files, for faster response
            if(toupper(input$clearData)=="CLEAR"){
                system(paste0("rm /pontus/data/shiny/emil/pbs",folderName,"*"))
            }
            
            ## preparing arguments for whole genome pbs calculations, or between and total for windows
            n<-length(chr)
            
            if(input$together=="YES"){
                wgPBSpos<-wholeGenomePBS(windows="NO",fall=fall,nInd=nInd,pos=pos,rs=rs,chr=chr,n=n,al12=numeric(n),al13=numeric(n),al23=numeric(n),bal12=numeric(n),bal13=numeric(n),bal23=numeric(n),pbs=numeric(n),pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,winSize=input$winSize,minWin=input$minWin,shinyPBS=shinyPBS,SNPsinChr=SNPsinChr[,1],PBSonly=T,maxChr=max(chr),FstOnly=input$FstOnly=="YES")
                wgPBSpos$pos<-pos
                wgPBSpos$chr<-chr
                winPBSpos<-extractWindow(wgPBSpos=wgPBSpos,start=start,end=end,chr=input$chr,minWin=input$minWin,ifWindows="NO",FstOnly=input$FstOnly=="YES")
                winPos<-winPBSpos[["winPos"]]
                winPBS<-winPBSpos[[ifelse(input$FstOnly=="YES","Fst12","winPBS")]]
                
                wgPBSpos2<-wholeGenomePBS(windows="YES",fall=fall,nInd=nInd,pos=pos,rs=rs,chr=chr,n=n,al12=numeric(n),al13=numeric(n),al23=numeric(n),bal12=numeric(n),bal13=numeric(n),bal23=numeric(n),pbs=numeric(n),pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,winSize=input$winSize,minWin=input$minWin,shinyPBS=shinyPBS,SNPsinChr=SNPsinChr[,1],maxChr=max(chr),FstOnly=input$FstOnly=="YES")
                winPBSpos2<-extractWindow(wgPBSpos=wgPBSpos2,start=start,end=end,chr=input$chr,minWin=input$minWin,ifWindows="YES",FstOnly=input$FstOnly=="YES")
                winPos2<-winPBSpos2[["winPos"]]
                winPBS2<-winPBSpos2[[ifelse(input$FstOnly=="YES","Fst12","winPBS")]]
            } else{
                
                wgPBSpos<-wholeGenomePBS(windows=input$ifWindows,fall=fall,nInd=nInd,pos=pos,rs=rs,chr=chr,n=n,al12=numeric(n),al13=numeric(n),al23=numeric(n),bal12=numeric(n),bal13=numeric(n),bal23=numeric(n),pbs=numeric(n),pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,winSize=input$winSize,minWin=input$minWin,shinyPBS=shinyPBS,SNPsinChr=SNPsinChr[,1],PBSonly=(input$ifWindows=="NO"),maxChr=max(chr),FstOnly=input$FstOnly=="YES")
              
                if(input$ifWindows=="NO"){
                  wgPBSpos$pos<-pos
                  wgPBSpos$chr<-chr
                }
                winPBSpos<-extractWindow(wgPBSpos=wgPBSpos,start=start,end=end,chr=input$chr,minWin=input$minWin,ifWindows="NO",FstOnly=input$FstOnly=="YES")
                winPos<-winPBSpos[["winPos"]]
                winPBS<-winPBSpos[[ifelse(input$FstOnly=="YES","Fst12","winPBS")]]
              
            }   
            ## function that extracts the chosen interval of the chosen chromosome

            
            if(length(winPos)==0&length(winPBS)==0){
                plot(0:1,0:1,main="no SNPs in interval selected!")
                return()
            }
            
            ## for displaying max function in chosen interval
          
            mmax<-max(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]],na.rm=T)
            cexaxis<-1.5
            cexlab<-2
            cexmain<-2
            cexpoint<-1.5
            
            ## plots PBS as function of position, for windows midpos of window
            ## put gene names in plot if 2 Mb or below
            if(max(winPos)-min(winPos)<=2.1e6){
                if(input$together=="YES"){
                    par(mfrow=c(1,3))
                    nf <- layout(matrix(c(1,2,3), 3, 1, byrow=TRUE),heights=c(1,0.5),widths=c(8))
                    par(mar=c(0.2, 6.2, 3, 6.3) + 0.1)
                    
                    plot(winPos/1e6,winPBS,type="p",pch=16,col="darkgrey",xaxt="n",ylab=ifelse(input$FstOnly=="YES","Fst12","PBS"),cex=cexpoint,cex.axis=cexaxis,cex.lab=cexlab,cex.main=cexmain)
                    abline(v=input$posPBS/1e6)
                    
                    par(new=TRUE)
                    plot(winPos2/1e6,winPBS2,type="l",pch=17,col="darkred",xaxt="n",yaxt="n",ylab="",cex=cexpoint,cex.axis=cexaxis,cex.lab=cexlab,cex.main=cexmain)
                    
                    ## legend has to emulate behavior of plot commands
                    ## perhaps find out if windows or not and plot accordingly so windows always red, lines and y-axis on the right
                    
                    axis(4,cex=cexpoint,cex.axis=cexaxis,cex.lab=cexlab,cex.main=cexmain)
                    mtext(paste0(ifelse(input$FstOnly=="YES","Fst12","PBS"))," red points",side=4,line=3,cex=1.3)
                    legend("topright",col=c("darkred","darkgrey"),lty=1,pch=c(16,17),legend=c("YES","NO"),cex=1.5)
                    ## will put x % quantile line in plot
                    if(input$quantileLine>0 & input$quantileLine<100) {
                        abline(h=quantile(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]],input$quantileLine/100,na.rm=T),col="red")  
                    }
                    ## function for putting genes in plot
                    plotGenes(input$chr,min(winPos)/1e6,max(winPos)/1e6,geneFile)
                    
                } else{
                    par(mfrow=c(1,3))
                    nf <- layout(matrix(c(1,2,3), 3, 1, byrow=TRUE),heights=c(1,0.5),widths=c(8))
                    par(mar=c(0.2, 6.2, 3, 6.3) + 0.1)
                    plot(winPos/1e6,winPBS,type=input$plotType,pch=16,col="darkblue",xaxt="n",ylab=ifelse(input$FstOnly=="YES","Fst12","PBS"),main=paste("max ",ifelse(input$FstOnly=="YES","Fst12","PBS")," value in interval is pos ",winPos[which(winPBS==max(winPBS,na.rm=T))]," with ",round(winPBS[which(winPBS==max(winPBS,na.rm=T))],digits = 3),sep=""),cex=cexpoint,cex.axis=cexaxis,cex.lab=cexlab,cex.main=cexmain)
                    abline(v=input$posPBS/1e6)
                    ## will put x % quantile line in plot
                    if(input$quantileLine>0 & input$quantileLine<100) {
                        abline(h=quantile(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]],input$quantileLine/100,na.rm=T),col="red")  
                    }
                    ## function for putting genes in plot
                    plotGenes(input$chr,min(winPos)/1e6,max(winPos)/1e6,geneFile)
                }
                
            } else{
                if(input$together=="YES"){
                    par(mfrow=2:1)
                    par(mar=c(5.1, 6.2, 3, 6.3))
                    
                    plot(winPos/1e6,winPBS,type="p",pch=16,col="darkgrey",xaxt="n",ylab=ifelse(input$FstOnly=="YES","Fst12","PBS"),cex=cexpoint,cex.axis=cexaxis,cex.lab=cexlab,cex.main=cexmain)
                    abline(v=input$posPBS/1e6)
                    
                    par(new=TRUE)
                    plot(winPos2/1e6,winPBS2,type="l",pch=17,col="darkred",xaxt="n",yaxt="n",ylab="",cex=cexpoint,cex.axis=cexaxis,cex.lab=cexlab,cex.main=cexmain)
                    
                    ## legend has to emulate behavior of plot commands
                    ## perhaps find out if windows or not and plot accordingly so windows always red, lines and y-axis on the right
                    
                    axis(4,cex=cexpoint,cex.axis=cexaxis,cex.lab=cexlab,cex.main=cexmain)
                    mtext(paste0(ifelse(input$FstOnly=="YES","Fst12","PBS"))," red points",side=4,line=3,cex=1.3)
                    legend("topright",col=c("darkred","darkgrey"),lty=1,pch=c(16,17),legend=c("YES","NO"),cex=1.5)
                    ## will put x % quantile line in plot
                    if(input$quantileLine>0 & input$quantileLine<100) {
                        abline(h=quantile(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]],input$quantileLine/100,na.rm=T),col="red")  
        }
                    
                }
                else{
                    par(mfrow=2:1)
                    par(mar=c(5.1, 6.2, 3, 6.3))
                    plot(winPos/1e6,winPBS,type=input$plotType,xlab=paste("Position (Mb) on ",input$chr,sep=""),pch=16,col="darkblue",ylab=ifelse(input$FstOnly=="YES","Fst12","PBS"),main=paste("max ",ifelse(input$FstOnly=="YES","Fst12","PBS")," value in interval is pos ",winPos[which(winPBS==max(winPBS,na.rm=T))]," with ",round(winPBS[which(winPBS==max(winPBS,na.rm=T))],digits = 3),sep=""),cex=cexpoint,cex.axis=cexaxis,cex.lab=cexlab,cex.main=cexmain)
                    abline(v=input$posPBS/1e6)
                    if(input$quantileLine>0 & input$quantileLine<100) {
                        abline(h=quantile(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]],input$quantileLine/100,na.rm=T),col="red")  
                    }
                }
            }
    
            ## genomewide PBS values histogram, with max PBS value and max PBS value in interval as vertical lines
            par(mar=c(5.1, 6.2, 3, 6.3))
            hist(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]],br=100,col="darkred",xlab=ifelse(input$FstOnly=="YES","Fst12","PBS"),main=paste0("distribution of ",ifelse(input$FstOnly=="YES","Fst12","PBS")," values whole genome (blue max in interval, red chosen quantile line)"),cex.axis=cexaxis,cex.lab=cexlab,cex.main=cexmain)
            ##abline(v=mmax,lwd=2,col="red")
            ##abline(v=max(winPBS),lwd=2,col="blue")
            
            abline(v=max(winPBS,na.rm=T),lwd=2,col="darkblue")
            abline(v=quantile(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]],input$quantileLine/100,na.rm=T),col="darkred",lwd=2)  
            legend("topright",c(paste("Top region\n(",round(mean(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]]<max(winPBS,na.rm=T),na.rm=T)*100,2),"% )",sep=""),paste(input$quantileLine,"%")),lwd=1,col=c("darkblue","darkred"),cex=cexlab)
            
            

        })


    output$PBSsinglePlot <- renderPlot({
      withProgress(PBSsinglePlot(),message="Calculation in progress, PROGESS BAR DOES NOT MEAN ANYTHING!!!!")
        

        })
 
 #########################################
 
 manhattanPlot <- eventReactive(input$runManhattan, {
          
     if(pw(input,passWord)){
       return()
     }
     if(input$pop1==input$pop2 | input$pop1==input$pop3 | input$pop2==input$pop3 ){
       return(plot(1,main="You have picked on population at least twice!"))
     }
     ##input <- list(start=57,end=-2,chr=4,ifWindows="YES",winSize=50000,pop1="Ghana",pop2="Kenya",pop3="Uganda",minWin=10,together="YES",quantileLine=95,FstOnly="YES")

     ## clearing saved shiny saved files, for faster response
     if(toupper(input$clearData)=="CLEAR"){
       system(paste0("rm /pontus/data/shiny/emil/pbs",folderName,"*"))
     }
     
     ##preparing arguments for whole genome pbs calculations, or between and total for windows
     n<-length(chr)
     
     ## calculating PBS values for whole genome, either windows or single marker, C++ functions for speed
     wgPBSpos<-wholeGenomePBS(windows=input$ifWindows,fall=fall,nInd=nInd,pos=pos,rs=rs,chr=chr,n=n,al12=numeric(n),al13=numeric(n),
                              al23=numeric(n),bal12=numeric(n),bal13=numeric(n),bal23=numeric(n),pbs=numeric(n),
                              pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,winSize=input$winSize,minWin=input$minWin,shinyPBS=shinyPBS,SNPsinChr=SNPsinChr[,1],PBSonly=(input$ifWindows=="NO"),maxChr=max(chr),FstOnly=(input$FstOnly=="YES"))
     
     ## remove those with PBS lower than 0.4
     
     ## plots PBS as function of position, for windows midpos of window
     ## put gene names in plot if 2 Mb or below
     
     manPBS<-wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]][which(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]]>0.1)]
     if(input$ifWindows=="YES"){
      manChr<-wgPBSpos$chr[which(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]]>0.1)]
     } else{
       manChr<-chr[which(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]]>0.1)]
     }
     
    ccol <- c("lightblue","darkblue")
    
    plot(manPBS,col=ccol[manChr%%2+1],xaxt="n",pch=16,ylab=ifelse(input$FstOnly=="YES","Fst12","PBS"),xlab="chromosomes")
    tabC <- table(manChr)
    mm <- cumsum(tabC)-tabC/2
    mtext(1:29,1,1:29%%3,at=unlist(mm))
    abline(h=quantile(wgPBSpos[[ifelse(input$FstOnly=="YES","Fst12","PBS")]],input$quantileLine/100,na.rm=T),col="red")  
    
       
     
   ##}) 
 })

    output$manhattanPlot <- renderPlot({
      withProgress(manhattanPlot(),message="Calculation in progress, PROGESS BAR DOES NOT MEAN ANYTHING!!!!")
        
        })
    

############################################


 runTable <- eventReactive(input$runTable, {
     ## setwd("/home/emil/shiny/waterbuckPBSV3test/")
     ##input <- list(start=63,end=65,chr=13,ifWindows="NO",winSize=50000,pop1="Ghana",pop2="Kenya",pop3="Uganda",minWin=10,FstOnly="NO")
     ##shinyPBS <- "~/tmp"
   
     if(pw(input,passWord,type="matrix")){
      return()
     }
   
    if(input$pop1==input$pop2 | input$pop1==input$pop3 | input$pop2==input$pop3 ){
      
      return(data.frame("You have picked on population at least twice!"))
    }
   
     start<-input$start*1e06
     end<-input$end*1e06
     
     if(toupper(input$clearData)=="CLEAR"){
       system(paste0("rm /pontus/data/shiny/emil/pbs",folderName,"*"))
     }
     
     n<-length(chr)
     
     wgTable<-wholeGenomePBS(windows=input$ifWindows,fall=fall,nInd=nInd,pos=pos,rs=rs,chr=chr,n=n,al12=numeric(n),al13=numeric(n),
                             al23=numeric(n),bal12=numeric(n),bal13=numeric(n),bal23=numeric(n),pbs=numeric(n),
                             pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,winSize=input$winSize,minWin=input$minWin,shinyPBS=shinyPBS,SNPsinChr=SNPsinChr[,1],maxChr=max(chr),FstOnly=(input$FstOnly=="YES"))
        
     
     if(!(any(start<=wgTable$pos & end>=wgTable$pos & input$chr==wgTable$chr))){
       inInterval<-data.frame(x="no SNPs in interval selected!")
         ##dT<- dTable(df)
         return(inInterval)
     }
    
     ## function for generating the table outputted, either windows or single marker, either top 1000 PBS or whole genome (all parameter)
     
     if(input$FstOnly=="YES"){
       inInterval<-PBSTable(wgTable=wgTable,chr=input$chr,pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,ifWindows=input$ifWindows,winSize=input$winSize,minWin=input$minWin,genes=dat,all="NO",start=start,end=end,thisChr=input$chr,maxChr=max(chr),FstOnly=(input$FstOnly=="YES"))
     } else{
      inInterval<-PBSTable(wgTable=wgTable,chr=input$chr,pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,ifWindows=input$ifWindows,winSize=input$winSize,minWin=input$minWin,genes=dat,all="NO",start=start,end=end,thisChr=input$chr,maxChr=max(chr))
     }
     ##dT<- dTable(inInterval,sPaginationType = 'full_numbers' ,iDisplayLength = 50)  
     
     inInterval<-inInterval[ !is.na(as.numeric(inInterval[,ifelse(input$FstOnly=="YES","Fst12","PBS")])),]
     return(inInterval) 
        
  })
    
    output$tableTop <-  renderTable( {
        withProgress(runTable(),message="Calculation in progress, PROGESS BAR DOES NOT MEAN ANYTHING!!!!")
        
        
      })
      
      
    
########################

runPBSgenes <- eventReactive(input$runPBSgenes, {
    if(pw(input,passWord,type="matrix")){
      return(data.frame("Correct password Needed!"))
    }
  
    if(input$pop1==input$pop2 | input$pop1==input$pop3 | input$pop2==input$pop3 ){
        return(data.frame("You have picked on population at least twice!"))
    }
    if(toupper(input$clearData)=="CLEAR"){
      system(paste0("rm /pontus/data/shiny/emil/pbs",folderName,"*"))
    } 
    
    load("data/geneInfo.Rdata")
    ##input <- list(start=60,end=62,chr=11,ifWindows="NO",winSize=50000,pop1="Ghana",pop2="Kenya",pop3="Tanzania",minWin=10,together="YES",quantileLine=95,FstOnly="YES")
    ##shinyPBS <- "~/tmp/"
    
    start<-input$start*1e06
    end<-input$end*1e06
    
    ## if chr negative wholeGenome results
    
    ## pos in genes
    ## SNPs in each gene
    ## call C++ Fst calculator - perhaps do new R function with this, gives back table with gene info
    
    nGene<-length(unique(SNPsInGene$gene))
    nPos<-length(SNPsInGene$pos)
    
    posGene<-SNPsInGene[,"pos"]
    chrGene<-SNPsInGene[,"chr"]
    
    fallGene<-SNPsInGene[,c(input$pop1,input$pop2,input$pop3)]
    
    whichGene<-SNPsInGene$whichGeneIndex-1
  
    geneTable<-genePBS(fallGene=fallGene,nInd=nInd,posGene=posGene,chrGene=chrGene,nPos=nPos,start=start,end=end,inputChr=input$chr,al12=numeric(nPos),
                       al13=numeric(nPos),al23=numeric(nPos),bal12=numeric(nPos),bal13=numeric(nPos),bal23=numeric(nPos),pbs=numeric(nPos),
                       pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,minWin=input$minWin,whichGene=whichGene,nGene=nGene,geneRange=geneRange,wholeGenome=(input$chr<0),shinyPBS=shinyPBS,FstOnly=input$FstOnly=="YES")
    
    geneTable<-as.data.frame(geneTable,stringsAsFactors = F)
    
    if(input$FstOnly=="YES"){
      geneTable[,c("Fst12","chr","SNPsinGene","Fstquantile")]<-apply(geneTable[,c("Fst12","chr","SNPsinGene","Fstquantile")],2,as.numeric)
      geneTable<-geneTable[ order(geneTable$Fst12,decreasing = T),]
    } else{
      geneTable[,c("PBS","Fst12","Fst13","Fst23","chr","SNPsinGene","PBSquantile")]<-apply(geneTable[,c("PBS","Fst12","Fst13","Fst23","chr","SNPsinGene","PBSquantile")],2,as.numeric)
      geneTable<-geneTable[ order(geneTable$PBS,decreasing = T),]
    }
    
    
    if(nrow(geneTable)==0){
      geneTable<-data.frame(x="no genes in interval selected!")
      return(geneTable)
    }
    
    ##function for generating the table outputted, either windows or single marker, either top 1000 PBS or whole genome (all parameter)
    ##dT<- dTable(geneTable,sPaginationType = 'full_numbers' ,iDisplayLength = 50)
    geneTable<-geneTable[ !is.na(geneTable[,ifelse(input$FstOnly=="YES","Fst12","PBS")]),]
    return(geneTable) 
  })



output$genePBS <- renderTable( {
  withProgress(runPBSgenes(),message="Calculation in progress, PROGESS BAR DOES NOT MEAN ANYTHING!!!!")
    
    })

#######################

tableTopWg <- eventReactive(input$runWG, {
    if(pw(input,passWord,type="matrix")){
        return(data.frame("Correct password Needed!"))

    }
    if(input$pop1==input$pop2 | input$pop1==input$pop3 | input$pop2==input$pop3 ){
      return(data.frame("You have picked on population at least twice!"))
    }

    ##input <- list(start=60,end=-2,chr=11,ifWindows="YES",winSize=50000,pop1="MAN",pop2="TEEN",pop3="YRI",minWin=10)
    ##input <- list(start=60,end=66,chr=11,ifWindows="NO",winSize=50000,pop1="NAT",pop2="CEU",pop3="CHB",minWin=10)
    ##shinyPBS <- "~/tmp"
    ##start<-input$start*1e06
    ##end<-input$end*1e06
    
    if(toupper(input$clearData)=="CLEAR"){
      system(paste0("rm /pontus/data/shiny/emil/pbs",folderName,"*"))
    }
    
    
    n<-length(chr)
    
    wgTable<-wholeGenomePBS(windows=input$ifWindows,fall=fall,nInd=nInd,pos=pos,rs=rs,chr=chr,n=n,al12=numeric(n),al13=numeric(n),
                            al23=numeric(n),bal12=numeric(n),bal13=numeric(n),bal23=numeric(n),pbs=numeric(n),
                            pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,winSize=input$winSize,minWin=input$minWin,shinyPBS=shinyPBS,SNPsinChr=SNPsinChr[,1],maxChr=max(chr),FstOnly=input$FstOnly=="YES")
    
    
    ##if(!(any(start<=wgTable$pos & end>=wgTable$pos & input$chr==wgTable$chr))){
    ##  df<-data.frame(x="No SNPs in interval selected!")
    ##  dT<- dTable(df)
    ##  return(dT)
    ##}
    
    if(input$FstOnly=="YES"){
      top100 <-  order(as.numeric(wgTable$Fst12),decreasing=T)[1:100]
    } else{
      top100 =  order(as.numeric(wgTable$PBS),decreasing=T)[1:100]
    }
    ##inInterval<-PBSTable(wgTable=wgTable,chr=input$chr,pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,ifWindows=input$ifWindows,winSize=input$winSize,minWin=input$minWin,genes=dat,all="NO",start=start,end=end,thisChr=input$chr)
    if(input$ifWindows=="YES"){
      
      top100res = cbind(wgTable[[ifelse(input$FstOnly=="YES","Fst12","PBS")]][top100],wgTable$pos[top100],wgTable$chr[top100],wgTable$SNPsinWin[top100])
      colnames(top100res)=c(ifelse(input$FstOnly=="YES","Fst12","PBS"),"Pos","Chr","SNPsinWin")
      df <- as.data.frame(top100res)
      df[,3]=as.character(df[,3])
    }else{
      top100res = cbind(wgTable[[ifelse(input$FstOnly=="YES","Fst12","PBS")]][top100],wgTable$rs[top100],wgTable$pos[top100],wgTable$chr[top100])
      colnames(top100res)=c(ifelse(input$FstOnly=="YES","Fst12","PBS"),"rs","Pos","Chr")
      df <- as.data.frame(top100res)
      df[,3]=as.character(df[,3])
    }
   
    
    ##dT<- dTable(inInterval,sPaginationType = 'full_numbers' ,iDisplayLength = 50)
    ##dT <- dTable(dF ,sPaginationType = 'full_numbers' ,iDisplayLength = 100)    
    ##write.csv(dF,file="/home/albrecht/public/albrecht/open/tmp/idastmp.csv")#/home/ida/web/shiny/PBSwithGIv2/outdata/top1000.csv")
  
    df<-df[ !is.na(df[,ifelse(input$FstOnly=="YES","Fst12","PBS")]),]
    
    return(df) 

})


output$tableTopWg <- renderTable({
  withProgress(tableTopWg(),message="Calculation in progress, PROGESS BAR DOES NOT MEAN ANYTHING!!!!")
    

    })
    
    
mistake <- eventReactive(input$runError, {
    if(pw(input,passWord,type="matrix")){
        return(data.frame("Correct password Needed!"))
    }
    
      file <- list.files("/var/log/shiny-server/",pattern="watebuckPBSV2",full.names = T)
            
            ## for pasting the most recent log file (tells of which error) to the shiny page
            scan(file,what="adad")
            
            paste("LOOKHERE           ",scan,"          LOOKHERE",sep="")
            
        })
        
        
    output$mistake <- renderText( {
      
      mistake()
      
    })
 
#################################################################
    #input <- list(start=60,end=62,chr=11,ifWindows="NO",winSize=50000,pop1="NAT",pop2="CEU",pop3="CHB",minWin=10)

    output$downloadData <- downloadHandler(
       

        filename = function() { paste('chr',input$chr,"_",input$start,input$end,input$pop1,input$pop2,input$pop3,"windows",input$ifWindows,ifelse(input$ifWindows=="NO","",input$winSize),'.txt',sep="") },
        content = function(file) {
            
          ##if(pw(input,passWord)){
          ##  return()
          ##}
          start<-input$start*1e06
          end<-input$end*1e06
          
       
          
          n<-length(chr)
          
          wgTable<-wholeGenomePBS(windows=input$ifWindows,fall=fall,nInd=nInd,pos=pos,rs=rs,chr=chr,n=n,al12=numeric(n),al13=numeric(n),
                                  al23=numeric(n),bal12=numeric(n),bal13=numeric(n),bal23=numeric(n),pbs=numeric(n),
                                  pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,winSize=input$winSize,minWin=input$minWin,shinyPBS=shinyPBS,SNPsinChr=SNPsinChr[,1],maxChr=max(chr),FstOnly=(input$FstOnly=="YES"))
          
          
          if(!(any(start<=wgTable$pos & end>=wgTable$pos & input$chr==wgTable$chr))){
            inInterval<-data.frame(x="no SNPs in interval selected!")
            ##dT<- dTable(df)
            return(inInterval)
          }
          
          ## function for generating the table outputted, either windows or single marker, either top 1000 PBS or whole genome (all parameter)
          
          if(input$FstOnly=="YES"){
            inInterval<-PBSTable(wgTable=wgTable,chr=input$chr,pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,ifWindows=input$ifWindows,winSize=input$winSize,minWin=input$minWin,genes=dat,all="NO",start=start,end=end,thisChr=input$chr,maxChr=max(chr),FstOnly=(input$FstOnly=="YES"))
          } else{
            inInterval<-PBSTable(wgTable=wgTable,chr=input$chr,pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,ifWindows=input$ifWindows,winSize=input$winSize,minWin=input$minWin,genes=dat,all="NO",start=start,end=end,thisChr=input$chr,maxChr=max(chr))
          }
          ##dT<- dTable(inInterval,sPaginationType = 'full_numbers' ,iDisplayLength = 50)  
          inInterval<-inInterval[ !is.na(inInterval[,ifelse(input$FstOnly=="YES","Fst12","PBS")]),]
          write.table(inInterval,file,sep="\t",qu=F,row=F)
            
        }
      
        ##})
            
        
    )
##################

        output$downloadgenePBS <- downloadHandler(
       

        filename = function() { paste('chr',input$chr,"_",input$start,input$end,input$pop1,input$pop2,input$pop3,"_minSNPs",input$minWin,'genePBS.txt',sep="") },
        content = function(file) {
            
         ## if(pw(input,passWord)){
         ##    return()
         ## }.
        if(input$pop1!=input$pop2 & input$pop1!=input$pop3 & input$pop2!=input$pop3 ){
          
            load("data/geneInfo.Rdata")
            start<-input$start*1e06
            end<-input$end*1e06
            
            ## if chr negative wholeGenome results
            
            ## pos in genes
            ## SNPs in each gene
            ## call C++ Fst calculator - perhaps do new R function with this, gives back table with gene info
            
            nGene<-length(unique(SNPsInGene$gene))
            nPos<-length(SNPsInGene$pos)
            
            posGene<-SNPsInGene[,"pos"]
            chrGene<-SNPsInGene[,"chr"]
            
            fallGene<-SNPsInGene[,c(input$pop1,input$pop2,input$pop3)]
            
            whichGene<-SNPsInGene$whichGeneIndex-1
            
            geneTable<-genePBS(fallGene=fallGene,nInd=nInd,posGene=posGene,chrGene=chrGene,nPos=nPos,start=start,end=end,inputChr=input$chr,al12=numeric(nPos),
                               al13=numeric(nPos),al23=numeric(nPos),bal12=numeric(nPos),bal13=numeric(nPos),bal23=numeric(nPos),pbs=numeric(nPos),
                               pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,minWin=input$minWin,whichGene=whichGene,nGene=nGene,geneRange=geneRange,wholeGenome=(input$chr<0),shinyPBS=shinyPBS,FstOnly=input$FstOnly=="YES")
            
            geneTable<-as.data.frame(geneTable,stringsAsFactors = F)
            
            if(input$FstOnly=="YES"){
              geneTable[,c("Fst12","chr","SNPsinGene","Fstquantile")]<-apply(geneTable[,c("Fst12","chr","SNPsinGene","Fstquantile")],2,as.numeric)
              geneTable<-geneTable[ order(geneTable$Fst12,decreasing = T),]
            } else{
              geneTable[,c("PBS","Fst12","Fst13","Fst23","chr","SNPsinGene","PBSquantile")]<-apply(geneTable[,c("PBS","Fst12","Fst13","Fst23","chr","SNPsinGene","PBSquantile")],2,as.numeric)
              geneTable<-geneTable[ order(geneTable$PBS,decreasing = T),]
            }
            
            
            if(nrow(geneTable)==0){
              geneTable<-data.frame(x="no genes in interval selected!")
              return(geneTable)
            }
            
            ##function for generating the table outputted, either windows or single marker, either top 1000 PBS or whole genome (all parameter)
            ##dT<- dTable(geneTable,sPaginationType = 'full_numbers' ,iDisplayLength = 50)
            geneTable<-geneTable[ !is.na(geneTable[,ifelse(input$FstOnly=="YES","Fst12","PBS")]),]
               
            ##function for generating the table outputted, either windows or single marker, either top 1000 PBS or whole genome (all parameter)            
            ##dT<- dTable(geneTable,sPaginationType = 'full_numbers' ,iDisplayLength = 50)               
            write.table(geneTable,file,sep="\t",qu=F,row=F)
          }
            
        }
    
    )
    


#######################3
    wholeGenome <- eventReactive(input$runDownload, {
        
            
        if(pw(input,passWord,type="matrix")){
            return(data.frame("Correct password Needed!"))

        }
            
        
        ## name the file with populations chosen, if windows and window size
        name<-paste(input$pop1,input$pop2,input$pop3,"windows",input$ifWindows,ifelse(input$ifWindows=="NO","",input$winSize),ifelse(input$FstOnly=="YES","FstOnly",""),'.gz',sep="")
        
        outFile <- paste("/home/albrecht/public/open/tmp/",name,sep="")
        
        ## if out file already exists then tell that and return that
        if(file.exists(outFile)){
            df<-data.frame(x=paste("Results already made\nget data here popgen.dk/albrecht/open/tmp/",basename(outFile),sep=""))
            return(df)
        }
        
        start<-input$start*1e06
        end<-input$end*1e06
        n<-length(chr)
        
        ## pass on arguments to scripts that will calculate this
        res <- list(windows=input$ifWindows,fall=fall,nInd=nInd,pos=pos,rs=rs,chr=chr,n=n,al12=numeric(n),al13=numeric(n),al23=numeric(n),bal12=numeric(n),bal13=numeric(n),bal23=numeric(n),pbs=numeric(n),pop1=input$pop1,pop2=input$pop2,pop3=input$pop3,winSize=input$winSize,minWin=input$minWin,shinyPBS=shinyPBS,SNPsinChr=SNPsinChr[,1],genes=dat,inputChr=input$chr,FstOnly=input$FstOnly)
      
        save(res,file="/home/albrecht/public/open/tmp/tmp.Rdata")
        
        ## has to be calculated by external scripts, too time consuming for shiny
        system("Rscript doStuff.R &")
        
        ## tell where to find outputted results
        
        
        ##plot(1,1,main=paste("get data here popgen.dk/albrecht/open/tmp/",input$pop1,input$pop2,input$pop3,"windows",input$ifWindows,ifelse(input$ifWindows=="NO","",input$winSize),".gz\n takes 10 min. go have coffee",sep=""))

        df<-data.frame(x=paste("Takes 10 min, get data here popgen.dk/albrecht/open/tmp/",basename(outFile),sep=""))
        ##dT<- dTable(df)
        return(df)
        

    })
    
    
    output$wholeGenome <- renderTable({
      withProgress(wholeGenome(),message="Calculation in progress, PROGESS BAR DOES NOT MEAN ANYTHING!!!!")
      
      
    })

           
            
  })
 


    
