
## read in C++ functions that calculate Reynolds Fst
## another that does weighted average of single locus estimators Fst
## lastly a function that calculates single PBS and between and total genomic variance

source("scripts/C_Fst_calculatorV2.R")
source("scripts/C_reynolds_Fst.R")

## for getting a faster rank function, data.table::frank
library(data.table)

##print(dirname(sys.frame(1)$ofile))

folderName<-scan("folderName",what="ll")[1]


errorPlot<-function(txt,cex=2,error="Error"){
    plot(0:1,0:1,col="transparent",xlab="",ylab="",axes=F)
    text(0.5,0.5,txt,cex=cex)
    mtext(error,1:4,col="red",cex=2)
}


pw<-function(x,passWord,type="plot"){
    tex<-"Enter correct password"
    if(type=="plot")
        errorPlot(tex)
    if(type=="matrix")
        matrix(tex,1,1)
    if(type=="text")
        tex
    return(!x$pw%in%passWord)
}



## this function plots then genes for a certain interval on a chromosome
## it uses a hg19 gene list and some trick with par & plot
plotGenes<-function(chr,min.pos,max.pos,geneFile){
  
  defaultPar<-par()
  
  xforrs = 0.03
  regsz = 1.2
  width=22 
  
  rsqplus = 0.045
  rightylabxplus=0.05
  xmargin = 0.005
  cexaxis=2
  cexlab=2
  adj = 0
  
  dat<-read.table(geneFile,as.is=T,head=T,comment.char="")
  xx2 = dat[dat[,"chrom"]==paste("chr",chr,sep="") & dat[,"cdsStart"]<max.pos*1e6 & dat[,"cdsEnd"] >min.pos*1e6,]
  
  start = xx2$txStart
  end   = xx2$txEnd
  nams  = xx2$name2
  cnts  = xx2$exonCount
  
  par(mar=c(5.2, 6.2, -0.1, 6.3) +0.1)
  
  plot(c(0,0),c(0,0),type="n",xlim=c(min.pos-adj,max.pos+adj),ylim=c(-0.8,0.1),xlab="",xaxs="i",yaxt='n',ylab="",main="",cex.lab=2.6,cex.axis=cexaxis,tck=-0.05)
  mtext(1,text=paste("Position on chromosome ",chr," (Mb)",sep=""),line=3,cex=1)
  
  ord <- order(start)
  start    <- start[ord]
  end      <- end[ord]
  exoncnts <- cnts[ord]
  nams     <- nams[ord]
  keep <- !duplicated(nams)
  start    <- start[keep]
  end      <- end[keep]
  exoncnts <- cnts[keep]
  nams     <- nams[keep]
  ord <- ord[keep]
  he       <- rep(c(0,-0.18,-0.36,-0.54,-0.72),100)[1:length(nams)]-0.05
  
  if(length(start)>0){
    segments(start/1e6, he, end/1e6, he)
    keep = !duplicated(nams)
    sapply(1:sum(keep),function(x){text((end[keep][x]+start[keep][x])/2e6,he[keep][x]+0.08,bquote(italic(.(nams[keep][x]))),cex=cexlab-0.6)})
    estart = as.numeric(unlist(sapply(xx2$exonStarts[ord],function(y){strsplit(y,",")[[1]]})))/1e6
    eend = as.numeric(unlist(sapply(xx2$exonEnds[ord],function(y){strsplit(y,",")[[1]]})))/1e6
    rect(estart,rep(he,xx2$exonCount[ord])-0.01,eend, rep(he,xx2$exonCount[ord])+0.01,col="black")
  }
  
  par(mar=defaultPar$mar)
}


genePBS<-function(fallGene,nInd,posGene,chrGene,nPos,start,end,inputChr,al12,al13,al23,bal12,bal13,bal23,pbs,pop1,pop2,pop3,minWin,whichGene,nGene,geneRange,whichFst,wholeGenome=F,shinyDir,FstOnly=F){
    ## fallGene=fallGene;nInd=nInd;posGene=posGene;chrGene=chrGene;nPos=nPos;start=start;end=end;inputChr=input$chr;al12=numeric(nPos);
    ## al13=numeric(nPos);al23=numeric(nPos);bal12=numeric(nPos);bal13=numeric(nPos);bal23=numeric(nPos);pbs=numeric(nPos);
    ## pop1=input$pop1;pop2=input$pop2;pop3=input$pop3;whichGene=whichGene;genes=nGene;
    ## geneRange=geneRange;wholeGenome=(input$chr<0);minwin=10;

    
    newName<-paste(paste0(shinyDir,"/",folderName,"/",".",pop1),pop2,pop3,minWin,"chr","ALL",start,end,"genePBS",ifelse(FstOnly,"Fst",""),sep=".")
    ## only load when doing windows, because single SNP PBS is not slow
    if(file.exists(newName)){
        load(newName)
        if(inputChr>0){
            geneFst3<-geneFst3[ as.numeric(geneFst3$chr)==inputChr & as.numeric(sapply(geneFst3$range,function(x) unlist(strsplit(x,"-"))[1])) > start & as.numeric(sapply(geneFst3$range,function(x) unlist(strsplit(x,"-"))[2])) < end,]
            
        }        
        return(geneFst3)
        
    }
    else{
        wholeGenome<-TRUE
        
        switch(whichFst,
               Reynolds={        
                   pbs2<-pbsCalculator_forCpp_ReynoldFst(fallGene[,pop1],fallGene[,pop2],fallGene[,pop3],nInd[pop1],nInd[pop2],nInd[pop3],posGene,chrGene,
                                                         n=nPos,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs)
               },
               Neis={        
                   pbs2<-pbsCalculator_forCpp_NeiFst(fallGene[,pop1],fallGene[,pop2],fallGene[,pop3],posGene,chrGene,
                                                     n=nPos,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs)
               },
               Hudsons={        
                   pbs2<-pbsCalculator_forCpp_HudsonFst(fallGene[,pop1],fallGene[,pop2],fallGene[,pop3],nInd[pop1],nInd[pop2],nInd[pop3],posGene,chrGene,
                                                        n=nPos,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs)
               },
               
               {
                   print("Not valid Fst estimator selected")
               }
               )
        
        geneFst<-list()
        for(e in c("al12","al13","al23")){
            
            ## preparing arguments for whole genome Fst window calculations
            ## calculating PBS values for windows via weighted average of single locus estimators Fst         
            
            cal2<-geneFstCPP2(vBetween=as.numeric(pbs2[[e]]),vTotal=as.numeric(pbs2[[paste("b",e,sep="")]]),nSNP=rep(nPos,nPos),Vcombined=numeric(nGene),whichGene=whichGene,
                              vBetweenGenes=numeric(nGene),vTotalGenes=numeric(nGene),genes=rep(nGene,nPos))
                        
            geneRes<-geneFst2(nSNP=nPos,vBetween=as.numeric(pbs2[[e]]),vTotal=as.numeric(pbs2[[paste("b",e,sep="")]]),whichGene=whichGene,genes=nGene)
            
            ## give gene name and calculate Fst
            ## give range of gene as argument as well            
            cal3<-lapply(cal2, function(x) x[1:nGene])
            
            ## extracting the arguments needed for further calculations from the function
            geneFst[[sub("al","Fst",e)]]<-cal3[["Vcombined"]]
            
        }
        geneFst[["chr"]]<-geneRange[,"chr"]
    
        geneFst[["range"]]<-geneRange[,"range"]
        geneFst[["geneNames"]]<-geneRange[,"name2"]
        geneFst[["SNPsinGene"]]<-geneRange[,"numberSNPs"]
        
        geneFst[["txStart"]]<-geneRange[,"txStart"]
        geneFst[["txEnd"]]<-geneRange[,"txEnd"]
        
        geneFst[["PBS"]]<-(-log(1-as.numeric(geneFst[["Fst12"]])) + -log(1-as.numeric(geneFst[["Fst13"]])) - -log(1-as.numeric(geneFst[["Fst23"]]))) / 2
        
        geneFst<-lapply(geneFst, function(x) x[which(as.numeric(geneFst[["SNPsinGene"]]) >= minWin)])
        
        if(wholeGenome){
            inIntervalIndex<-!logical(length(geneFst[[1]]))
        } else{
            inIntervalIndex<-(as.numeric(geneFst[["txEnd"]]) > start & as.numeric(geneFst[["txStart"]]) < end & as.numeric(geneFst[["chr"]]) == inputChr)
        }
        if(FstOnly){
            geneFst[["Fstquantile"]] <- (data.table::frank(c(as.numeric(geneFst[["Fst12"]])))/length(geneFst[[1]]))*100
        } else{
            geneFst[["PBSquantile"]] <- (data.table::frank(c(as.numeric(geneFst[["PBS"]])))/length(geneFst[[1]]))*100
        }
        if(wholeGenome){
            geneFst3<-do.call(cbind, geneFst)
        } else{
            geneFst2<-lapply(geneFst, function(x) x[which(as.numeric(geneFst[["txEnd"]]) > start & as.numeric(geneFst[["txStart"]]) < end & as.numeric(geneFst[["chr"]]) == inputChr)]) 
            geneFst3<-do.call(cbind, geneFst2)
        }

        geneFst3<-geneFst3[ order(as.numeric(geneFst3[,"PBS"]),decreasing = T),]
        geneFst3<-cbind(PBS=geneFst3[,"PBS"],geneFst3[,c("Fst12","Fst13","Fst23","chr","range","geneNames","SNPsinGene",ifelse(FstOnly,"Fstquantile","PBSquantile"))])               

        geneFst3<-as.data.frame(geneFst3,stringsAsFactors=F)
        save(geneFst3,file=newName)
        if(inputChr>0){
            
            geneFst3<-geneFst3[ as.numeric(geneFst3$chr)==inputChr & as.numeric(sapply(geneFst3$range,function(x) unlist(strsplit(x,"-"))[1])) > start & as.numeric(sapply(geneFst3$range,function(x) unlist(strsplit(x,"-"))[2])) < end,]

        }     
        return(geneFst3) 
    }
}

## this one calculate PBS values for whole genome, either single marker or window,
## it uses the aforementioned C++ functions
## if windows it looks if file already exists and it loads if so
wholeGenomePBS<-function(windows,fall,nInd,pos,rs,chr,n,al12,al13,al23,bal12,bal13,bal23,pbs,pop1,pop2,pop3,winSize,minWin,shinyDir,SNPsinChr,whichFst,PBSonly=F,maxChr,FstOnly=F){

    if(windows=="NO"){
        ## calculating PBS values for whole genome
        
        switch(whichFst,
               Reynolds={        
                   pbs2<-pbsCalculator_forCpp_ReynoldFst(fall[,pop1],fall[,pop2],fall[,pop3],nInd[pop1],nInd[pop2],nInd[pop3],pos,chr,
                                                         n=n,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs)
               },
               Neis={        
                   pbs2<-pbsCalculator_forCpp_NeiFst(fall[,pop1],fall[,pop2],fall[,pop3],pos,chr,
                                                     n=n,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs)
               },
               Hudsons={
      	           pbs2<-pbsCalculator_forCpp_HudsonFst(fall[,pop1],fall[,pop2],fall[,pop3],nInd[pop1],nInd[pop2],nInd[pop3],pos,chr,
                                                      n=n,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs)

               },
               {
                   print("Not valid Fst estimator selected")
               }
               )
        
        if(FstOnly & PBSonly){
            Fst12<-ifelse(pbs2[["al12"]]/pbs2[["bal12"]]>0.99,0.99,pbs2[["al12"]]/pbs2[["bal12"]])
            return(list(Fst12=Fst12))
        } else if(PBSonly){
            return(list(PBS=pbs2[["pbs"]]))
        } 
                
        Fst12<-ifelse(pbs2[["al12"]]/pbs2[["bal12"]]>0.99,0.99,pbs2[["al12"]]/pbs2[["bal12"]])
        Fst13<-ifelse(pbs2[["al13"]]/pbs2[["bal13"]]>0.99,0.99,pbs2[["al13"]]/pbs2[["bal13"]])
        Fst23<-ifelse(pbs2[["al23"]]/pbs2[["bal23"]]>0.99,0.99,pbs2[["al23"]]/pbs2[["bal23"]])

        if(FstOnly){
            wgTable <- list(chr=pbs2[["chr"]],pos=pbs2[["pos"]],rs=rs,Fst12=Fst12,fall[,pop1],fall[,pop2],fall[,pop3])
            names(wgTable)<-c("chr","pos","rs","Fst12",pop1,pop2,pop3)
        } else{
            wgTable <- list(chr=pbs2[["chr"]],pos=pbs2[["pos"]],rs=rs,PBS=pbs2[["pbs"]],Fst12,Fst13,
                            Fst23,fall[,pop1],fall[,pop2],fall[,pop3])
            names(wgTable)<-c("chr","pos","rs","PBS",paste("Fst",pop1,pop2,sep=""),paste("Fst",pop1,pop3,sep=""),paste("Fst",pop2,pop3,sep=""),pop1,pop2,pop3)
        }
        

        return(wgTable)
        
    } else{
        ## calculating between and total genomic variance for whole genome (functions gives those and PBS values)

        switch(whichFst,
               Reynolds={        
                   pbs2<-pbsCalculator_forCpp_ReynoldFst(fall[,pop1],fall[,pop2],fall[,pop3],nInd[pop1],nInd[pop2],nInd[pop3],pos,chr,
                                                         n=n,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs)
               },
               Neis={        
                   pbs2<-pbsCalculator_forCpp_NeiFst(fall[,pop1],fall[,pop2],fall[,pop3],pos,chr,
                                                     n=n,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs)
               },
               Hudsons={
               pbs2<-pbsCalculator_forCpp_HudsonFst(fall[,pop1],fall[,pop2],fall[,pop3],nInd[pop1],nInd[pop2],nInd[pop3],pos,chr,		
                                                         n=n,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs)		
               },
               {
                   print("Not valid Fst estimator selected")
               }
               )
        
        
        winFst<-list()

        ## generate a file name for the settings to save
        newName<-paste(paste0(shinyDir,"/",folderName,"/",".",pop1),pop2,pop3,winSize,ifelse(FstOnly,"Fst",""),"table",sep=".")
        ## only load when doing windows, because single SNP PBS is not slow
        if(file.exists(newName)){
            load(newName)
            if(minWin!=10){
                ## had to to put minimum in window here, otherwise does not work
                wgTable2<-lapply(wgTable,function(x) x[which(wgTable$SNPsinWin>minWin)])  
                return(wgTable2)
            } 
            return(wgTable)
            
        }
        else{
            
    for(e in c("al12","al13","al23")){
        
        ## preparing arguments for whole genome Fst window calculations
        nSNP<-rep(n,times=n)
        Vcombined<-numeric(n)
        regionSize<-integer(n)
        Vmidpos<-numeric(n)
        vBetween<-as.numeric(pbs2[[e]])   
        vTotal<-as.numeric(pbs2[[paste("b",e,sep="")]])      
        winSizeVector<-rep(winSize,times=n)
        SNPsinChrVector<-integer(n)
        SNPsinChrVector[1:maxChr]<-SNPsinChr
        maxChr<-as.numeric(max(chr))
        
        ## calculating PBS values for windows via weighted average of single locus estimators Fst
        cal2<-calculate_genomeFst_forCpp(pos,chr,vBetween,vTotal,nSNP=nSNP,Vcombined=Vcombined,Vmidpos=Vmidpos,winSize=winSizeVector,regionSize=regionSize,SNPsinChr=SNPsinChrVector,maxChr=maxChr)
        
        ## extracting the arguments needed for further calculations from the function
        winFst[[e]]<-cal2[["Vcombined"]]
        winFst[[paste("midPos_",e,sep="")]]<-cal2[["Vmidpos"]]
        winFst[[paste("regSize_",e,sep="")]]<-cal2[["regionSize"]]
        winFst[["chr"]]<-cal2[["chr"]]
        
    }
            
            winMatrix<-do.call(cbind, winFst)
            
            ## only windows with actual SNPs in them, because C++ inline functions return whole input argument
            winMatrix2<-winMatrix[ winMatrix[,"regSize_al12"]!=0&winMatrix[,"regSize_al13"]!=0&winMatrix[,"regSize_al23"]!=0,]
            
            ## only take those windows with SNPs and above a chosen minimal window size
            winMatrix2<-winMatrix2[ winMatrix2[,"regSize_al12"]>=minWin,]

            if(FstOnly){
                return( list(chr=winMatrix2[,"chr"],pos=winMatrix2[,"midPos_al12"],Fst12=winMatrix2[,"al12"],SNPsinWin=winMatrix2[,"regSize_al12"]))
            }
            
            tempPBS<-(-log(1-winMatrix2[,"al12"]) + -log(1-winMatrix2[,"al13"]) - -log(1-winMatrix2[,"al23"])) / 2
            wgTable <- list(chr=winMatrix2[,"chr"],pos=winMatrix2[,"midPos_al12"],PBS=tempPBS,SNPsinWin=winMatrix2[,"regSize_al12"],winMatrix2[,"al12"],winMatrix2[,"al13"],winMatrix2[,"al23"])
            
            names(wgTable)<-c("chr","pos","PBS","SNPsinWin",paste("Fst",pop1,pop2,sep=""),paste("Fst",pop1,pop3,sep=""),paste("Fst",pop2,pop3,sep=""))

            ## save whole genome results if did not exist already - for faster respone from shiny website
            save(wgTable,file=newName)
            return(wgTable)
        }
        
    }
}


## this function extracts the desired interval of the PBS & pos values
## if start is negative it means whole chromosome
## if end is negative it means end size interval around start
extractWindow<-function(wgPBSpos,start,end,chr,minWin,ifWindows,FstOnly=F){
    
    if(start<0){
        winPBS<-wgPBSpos[[ifelse(FstOnly,"Fst12","PBS")]][which(wgPBSpos[["chr"]]==chr)]
        winPos<-wgPBSpos[["pos"]][which(wgPBSpos[["chr"]]==chr)]
    }else if(end<0){
        
        winPBS<-wgPBSpos[[ifelse(FstOnly,"Fst12","PBS")]][which(wgPBSpos[["pos"]]>=start+end&wgPBSpos[["pos"]]<=start-end&wgPBSpos[["chr"]]==chr)]
        winPos<-wgPBSpos[["pos"]][which(wgPBSpos[["pos"]]>=start+end&wgPBSpos[["pos"]]<=start-end&wgPBSpos[["chr"]]==chr)]    
        
    }else{
        winPBS<-wgPBSpos[[ifelse(FstOnly,"Fst12","PBS")]][which(wgPBSpos[["pos"]]>=start&wgPBSpos[["pos"]]<=end&wgPBSpos[["chr"]]==chr)]
        winPos<-wgPBSpos[["pos"]][which(wgPBSpos[["pos"]]>=start&wgPBSpos[["pos"]]<=end&wgPBSpos[["chr"]]==chr)]    
        
    }
    if(FstOnly){
        return(list(Fst12=winPBS,winPos=winPos))
    } else{
        return(list(winPBS=winPBS,winPos=winPos))
    }
}



## calculate closests gene to single marker
closestGene<-function(chr,pos,genes){
    genes2<-genes[ genes$chr==chr,]
    ## find which distance is smallest, looks both at cdsStart & cdsEnd
    minDist<-sort.list(c(abs(pos-genes2[ ,"cdsStart"]),abs(pos-genes2[ ,"cdsEnd"])))[1] 
    ## if closer to cdsEnd, will be in second half of list: and will have index of original index + nrow(genes)
    ## because it is sorted as minDist of both cdsStart and cdsEnd distances, so twice as many as genes
    minDist2<-ifelse(minDist>nrow(genes2),minDist-nrow(genes2),minDist) 
    return(genes2[minDist2 ,"name2"])
}

## this functions does the PBS table with or without windodws,
## it puts Fst, PBS values, genes, quantile int it
PBSTable<-function(wgTable,chr,pop1,pop2,pop3,ifWindows,winSize,minWin,genes,all,start,end,thisChr,maxChr,FstOnly=F){
    ##wgTable=wgTable;chr=input$chr;pop1=input$pop1;pop2=input$pop2;pop3=input$pop3;ifWindows=input$ifWindows;winSize=input$winSize;genes=dat;all="NO";input$start=start;input$end=end
    wgTable$chr<-as.numeric(wgTable$chr)
    wgTable$pos<-as.numeric(wgTable$pos)
    if(FstOnly){
        wgTable$Fst12<-as.numeric(wgTable$Fst12)
    } else{
        wgTable$PBS<-as.numeric(wgTable$PBS)
    }
    
    ## to extract the right interval of values again end<0 interval around start, start<0 whole chromosome
    if(all=="YES"){
        inIntervalIndex<-rep(TRUE,length(wgTable$pos))
        notInIntervalIndex<-!(inIntervalIndex)
    } else if(start<0){
        inIntervalIndex<-which(wgTable$chr==chr)
    notInIntervalIndex<-which(wgTable$chr!=chr)
    } else if(end<0){
        inIntervalIndex<-which(wgTable$pos>=(start-abs(end))&wgTable$pos<=(start+abs(end))&wgTable$chr==chr)
        notInIntervalIndex<-which(wgTable$pos<(start-abs(end))|wgTable$pos>(start+abs(end))|wgTable$chr!=chr)    
    } else{
        inIntervalIndex<-which(wgTable$pos>=start&wgTable$pos<=end&wgTable$chr==chr)
        notInIntervalIndex<-which(wgTable$pos<start|wgTable$pos>end|wgTable$chr!=chr)
        
    }
        
    ## calculate which quantile PBS values are in, by taking their 1-(rank/nSNP), only for those in desired window
    PBSquantile <- data.table::frank(c(wgTable[[ifelse(FstOnly,"Fst12","PBS")]][inIntervalIndex],wgTable[[ifelse(FstOnly,"Fst12","PBS")]][notInIntervalIndex]))[1:length(inIntervalIndex)]/length(wgTable[[1]])
    
    ## different things have to be displayed depending on if with windows or not
    if(ifWindows=="NO"){
        if(FstOnly){
            inInterval<-cbind(wgTable[["chr"]][inIntervalIndex],wgTable[["pos"]][inIntervalIndex],wgTable[["rs"]][inIntervalIndex],wgTable[["Fst12"]][inIntervalIndex],wgTable[[pop1]][inIntervalIndex],wgTable[[pop2]][inIntervalIndex],PBSquantile*100)
            colnames(inInterval)<-c(names(wgTable)[1:6],"quantile")
            inInterval[,"Fst12"]<-as.numeric(inInterval[,"Fst12"])
        } else{
            inInterval<-cbind(wgTable[["chr"]][inIntervalIndex],wgTable[["pos"]][inIntervalIndex],wgTable[["rs"]][inIntervalIndex],wgTable[["PBS"]][inIntervalIndex],wgTable[[paste("Fst",pop1,pop2,sep="")]][inIntervalIndex],wgTable[[paste("Fst",pop1,pop3,sep="")]][inIntervalIndex],wgTable[[paste("Fst",pop2,pop3,sep="")]][inIntervalIndex],wgTable[[pop1]][inIntervalIndex],wgTable[[pop2]][inIntervalIndex],wgTable[[pop3]][inIntervalIndex],PBSquantile*100)
            colnames(inInterval)<-c(names(wgTable),"quantile")
            inInterval[,paste("Fst",pop1,pop2,sep="")]<-as.numeric(inInterval[,paste("Fst",pop1,pop2,sep="")])
            inInterval[,paste("Fst",pop1,pop3,sep="")]<-as.numeric(inInterval[,paste("Fst",pop1,pop3,sep="")])
            inInterval[,paste("Fst",pop2,pop3,sep="")]<-as.numeric(inInterval[,paste("Fst",pop2,pop3,sep="")])
            inInterval[,"PBS"]<-as.numeric(inInterval[,"PBS"])
        }
        
    } else{
        
        if(FstOnly){
            inInterval<-cbind(wgTable[["chr"]][inIntervalIndex],wgTable[["pos"]][inIntervalIndex],wgTable[["Fst12"]][inIntervalIndex],wgTable[["SNPsinWin"]][inIntervalIndex],PBSquantile*100)
            colnames(inInterval)<-c(names(wgTable),"quantile")
            inInterval[,"Fst12"]<-as.numeric(inInterval[,"Fst12"])

        } else{
            inInterval<-cbind(wgTable[["chr"]][inIntervalIndex],wgTable[["pos"]][inIntervalIndex],wgTable[["PBS"]][inIntervalIndex],wgTable[["SNPsinWin"]][inIntervalIndex],wgTable[[paste("Fst",pop1,pop2,sep="")]][inIntervalIndex],wgTable[[paste("Fst",pop1,pop3,sep="")]][inIntervalIndex],wgTable[[paste("Fst",pop2,pop3,sep="")]][inIntervalIndex],PBSquantile*100)
            colnames(inInterval)<-c(names(wgTable),"quantile")
            inInterval[,paste("Fst",pop1,pop2,sep="")]<-as.numeric(inInterval[,paste("Fst",pop1,pop2,sep="")])
            inInterval[,paste("Fst",pop1,pop3,sep="")]<-as.numeric(inInterval[,paste("Fst",pop1,pop3,sep="")])
            inInterval[,paste("Fst",pop2,pop3,sep="")]<-as.numeric(inInterval[,paste("Fst",pop2,pop3,sep="")])
            inInterval[,"PBS"]<-as.numeric(inInterval[,"PBS"])
        }
        
    }
    inInterval<-as.data.frame(inInterval,stringsAsFactors = F)

    keep<-which(!is.na(as.numeric(inInterval[,ifelse(FstOnly,"Fst12","PBS")])))
    inInterval<-inInterval[ keep,]
    
    inInterval[,"chr"]<-as.numeric(inInterval[,"chr"])
    inInterval[,"pos"]<-as.numeric(inInterval[,"pos"])
    inInterval[,"quantile"]<-as.numeric(inInterval[,"quantile"])
    
    ## if not whole genome only returns top 1000 PBS
    if(all!="YES"){
        inInterval<-inInterval[order(as.numeric(inInterval[,ifelse(FstOnly,"Fst12","PBS")]),decreasing=T),]
        inInterval<-inInterval[1:min(1000,nrow(inInterval)),]        
    }
    
    genes2<-lapply(1:maxChr,function(x) genes[genes$chr==x,])
    genes2<-lapply(1:maxChr,function(x) genes2[[x]][order(genes2[[x]]$cdsStart),] )
    names(genes2)<-1:maxChr
    
    if(ifWindows=="NO"){
        ## put AF in allele frequency column
        if(FstOnly){
            colnames(inInterval)[c(ncol(inInterval)-2,ncol(inInterval)-1)]<-sapply(colnames(inInterval)[c(ncol(inInterval)-2,ncol(inInterval)-1)],function(x) paste("AF_",x,sep=""))
        }
        else{
            colnames(inInterval)[c(ncol(inInterval)-3,ncol(inInterval)-2,ncol(inInterval)-1)]<-sapply(colnames(inInterval)[c(ncol(inInterval)-3,ncol(inInterval)-2,ncol(inInterval)-1)],function(x) paste("AF_",x,sep=""))
        }
        if(all=="YES"){
            ## precaculated closest gene for each marker
            load("data/closestGeneSingleMarker.Rdata") # object called closestGeneSingleMarker
            closestGeneSingleMarker<-closestGeneSingleMarker[keep,]
            inInterval$closestsGene<-closestGeneSingleMarker[,"gene"]
            colnames(inInterval)[ncol(inInterval)]<-"closests_gene"
        } else{
            inInterval$closestsGene<-apply(inInterval,1,function(x) closestGene(as.numeric(x[1]),as.numeric(x[2]),genes))
            colnames(inInterval)[ncol(inInterval)]<-"closests_gene"
        }
        
    } else{
        ## pos is actually midPos, so edges midPos +- winSize/2
        inInterval<-inInterval[ inInterval$SNPsinWin>minWin,]

        ## this function calculates intersecting genes of windows
        superFun2 <- function(thisChr){          
            thisPos <- as.numeric(inInterval[inInterval[,"chr"]==thisChr,"pos"])            
            ## for each chr it keeps track of which genes added, dynamic programming FTW
            minGene<-1
            maxGene<-1
            maxMaxGene<-length(genes2[[thisChr]]$chr)
            
            res<-character(length(thisPos))
            for(i in 1:length(thisPos)){                
                added<-F  
                res2<-""
                exclude<-c()
                sel<-c()
                ## this counts up minGene and maxGene for all genes left/downstream of window
                while(genes2[[thisChr]]$cdsEnd[minGene] < (thisPos[i]-winSize/2) &  maxGene < nrow(genes2[[thisChr]])){
                    ## minGene is first gene intersected by window
                    minGene <- minGene + 1
                    ## maxGene is first gene intersected by window
                    maxGene <- maxGene + 1
                    if(maxGene >= maxMaxGene){
                        maxGene<-maxMaxGene
                        break
                    }
                }          
                ## this counts up maxGene for genes left/downstream of window, genes intersecting window with cdsEnd, genes overlapping windows completely
                while((genes2[[thisChr]]$cdsEnd[maxGene] < (thisPos[i]-winSize/2) &  maxGene < nrow(genes2[[thisChr]])) 
                      | (genes2[[thisChr]]$cdsStart[maxGene] < (thisPos[i]-winSize/2) & genes2[[thisChr]]$cdsEnd[maxGene] > (thisPos[i]-winSize/2)) 
                      | (genes2[[thisChr]]$cdsStart[maxGene] < (thisPos[i]-winSize/2) & genes2[[thisChr]]$cdsEnd[maxGene] > (thisPos[i]+winSize/2))
                      &  maxGene < nrow(genes2[[thisChr]])){
                          
                          ## looks at if windows cdsEnd is left of window, means should not be added
                          if(genes2[[thisChr]]$cdsEnd[maxGene] < (thisPos[i]-winSize/2)){
                              exclude<-c(exclude,maxGene)
              maxGene<-maxGene+1
            } else{
              maxGene <- maxGene + 1
              if(maxGene >= maxMaxGene){
                maxGene<-maxMaxGene
                break
              }
              added<-T
            }
                      }                
                ## counts up maxGene if gene is in window or hads cdsStart intersecting window      
                while((genes2[[thisChr]]$cdsStart[maxGene] > (thisPos[i]-winSize/2) & genes2[[thisChr]]$cdsEnd[maxGene] < (thisPos[i]+winSize/2)) |
                      (genes2[[thisChr]]$cdsStart[maxGene] < (thisPos[i]+winSize/2) & genes2[[thisChr]]$cdsEnd[maxGene] > (thisPos[i]+winSize/2))
                & maxGene < nrow(genes2[[thisChr]])){
            
            maxGene <- maxGene + 1
            if(maxGene >= maxMaxGene){
                maxGene<-maxMaxGene
                break
            }
                    added<-T
                }
                
                ## looks if any genes intersect window
                if(added){
                    
                    ## because maxGene is one too far to the right of window
                    maxGene<-maxGene-1
                    sel<-minGene:maxGene
              
              ## checks if exclude is null
              if(is.null(exclude)){
                  res2<- unique(genes2[[thisChr]][sel,"name2"])
              } else{
                res2<- unique(genes2[[thisChr]][sel[!(sel%in%exclude)],"name2"])  
              }
              
              ## goes back to first gene not added, so none are missed for new window
              if(minGene>=2){
                  minGene<-minGene-1
                  maxGene<-minGene  
              } else{
                  ##this is if first gene found is added and then minGene should not go back to 0
                  maxGene<-minGene  
            }
              
          }
                
                res[i]<-paste(res2,collapse="/")
            }
            return(res)
      }
       
        if(all=="YES"){
            ## if whole genome does multithreaded apply for speed!
            if(require(parallel,character.only = TRUE)){
                all2 <- parallel::mclapply(1:maxChr,superFun2,mc.cores=maxChr)
                inInterval$genes_contained_in_Window <- unlist(all2)
            } else{
                print("calculating whole genome NOT threaded - might take long")
                all2 <- lapply(1:maxChr,superFun2)
                inInterval$genes_contained_in_Window <- unlist(all2)
            }
        } else{
            inInterval<-inInterval[ order(inInterval$pos),]
            all2<-sapply(chr,superFun2)       
            inInterval$genes_intersected_in_Window<-unlist(all2)
            inInterval<-inInterval[ order(as.numeric(inInterval[,ifelse(FstOnly,"Fst12","PBS")]),decreasing = T),]
        }
        
    }
    inInterval$pos<-as.character(inInterval$pos)    
    return(inInterval)
}

