## du -cksh *

########### do not change ################3
l<-commandArgs(TRUE)
getArgs<-function(x,l){
    unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]
}
Args<-function(l,args){
    if(length(l)%%2!=0){
        cat("Error -> not all options have an argument!\n")
        q("no")       
    }
    if(sum(grepl("^-",l))!=length(l)/2){
        cat("Error -> not correct number of options - or argument starting with -\n")
        q("no")       
    }
    ## this assumes -option argument structure!
    l<-sapply(seq(1,(length(l)-1),2),function(x) paste0(gsub("^-","",l[x]),"=",l[x+1]))
    
    if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
        cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument")
        q("no")
    }
    arguments<-list()
    for(a in names(args)){
        arguments[[a]]<-getArgs(a,l)
    }
    
    ## check for freqFolder or Qfile + Ffile file
    ##if("freqFolder"%in%names(arguments) & !all(c("famFile","bimFile")%in%names(arguments))){
    ##   cat("Error -> either freqFolder or famFile, bimFile or geneFile argument(s) is/are missing!\n")
    ##    q("no")  
    ##} else if("filterFile"%in%names(arguments) & !all(c("Qfile","Ffile","popsFile")%in%names(arguments))){
    ##   cat("Error -> either filterFile or Qfile, Ffile or popsFile argument(s) is/are missing!\n")
    ##    q("no")  
    ##} else
    ## check if NGSadmix or ADMIXTURE or plink
    if((!all(c("Qfile","Ffile","popsFile","filterFile")%in%names(arguments)) | !all(c("Qfile","Ffile","popsFile","bimFile")%in%names(arguments))) & !all(c("freqFolder","famFile","bimFile")%in%names(arguments))){
        cat("Error -> Qfile, Ffile, popsFile or freqFolder argument(s) is/are missing!\n")
        q("no")
    ## checks that not both NGSadmix and ADMIXTURE    
    } else if(all(c("Qfile","Ffile","popsFile","filterFile","bimFile")%in%names(arguments))){
        cat("Error -> both filterFile and bimFile have been specified!\n")
        q("no")         
    } else if(!c("folderName")%in%names(arguments)){
        cat("Error -> folderName has to be specified!\n")
        q("no")  
    } else if(!c("shinyDir")%in%names(arguments)){
        cat("Error -> shinyDir has to be specified!\n")
        q("no")    
    } else if(!c("fileGenes")%in%names(arguments)){
        cat("Error -> fileGenes has to be specified!\n")
        q("no")
    }
    if(any(!names(args)%in%names(arguments)&sapply(args,is.null))){
        cat("Error -> ",names(args)[!names(args)%in%names(arguments)&sapply(args,is.null)]," is not optional!\n")
        q("no")
    }    
    for(a in names(args))
        if(is.null(arguments[[a]]))
            arguments[[a]]<-args[[match(a,names(args))]]
    return(arguments)
}

print.args<-function(args,des){
    if(missing(des)){
        des<-as.list(rep("",length(args)))
        names(des)<-names(args)
    }
    cat("->  needed arguments:\n")
    mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
    cat("->  optional arguments (defaults):\n")
    mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
    q("no")
}
###### ####### ###### ###### ###### #######
## choose your parameters and defaults
## NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments

## getting arguments for run
    
args<-list(freqFolder="",
           famFile="",
           bimFile="",
           Qfile = "",
           Ffile = "",
           filterFile = "",
           popsFile= "",
           fileGenes=NULL,
           folderName=NULL,
           shinyDir=NULL               
           )
## if no argument aree given prints the need arguments and the optional ones with default
des<-list(freqFolder="folder with .frq files of the pops to be included - all .freq files in folder will be read",
          famFile=".fam file of the data analysed with --freq",
          bimFile=".bim file of the data analysed with --freq or ADMIXTURE",
          Qfile = "Q file from either ADMIXTURE or NGSadmix",
          Ffile = "F file from either ADMIXTURE or NGSadmix",
          filterFile = ".filter file, which is list of SNPs included in the analysis with NGSadmix",
          popsFile= "list of populations to be analysed, order should be same as columns in ADMIXTURE or NGSadmix",
          fileGenes="file with a list of genes and their range",
          folderName="name of Shiny and folder",
          shinyDir="name of folder to store data"              
          )

######################################
#######get arguments and add to workspace
### do not change
if(length(l)==0){
    print.args(args,des)
}
attach(Args(l,args))
args <- commandArgs(TRUE)
if(length(args)==0){
    cat("Arguments: output prefix\n")
    q("no")
}
###################################


    
## Have to be in folder when running!
setwd(".")


system("mkdir -p data",intern=T)
system(paste0("mkdir -p ",shinyDir))
system(paste0("mkdir -p ",shinyDir,"/",folderName),intern=T)

write(folderName,"folderName")
write(shinyDir,"shinyDir")

if(freqFolder!=""){
    
    freqFiles<-list.files(freqFolder,full=T,pattern=".frq$")
    
    l<-list()
    for(file in freqFiles){
        
        name<-unlist(strsplit(tail(unlist(strsplit(file,"/")),1),".frq"))[1]
        l[[name]]<-read.table(file,h=T,as.is=T)
        print("done:")
        print(file)
        
    }
    
    bim<-read.table(bimFile,as.is=T,colC=c("integer","character")[c(1,2,1,1,2,2)])
    fam<-read.table(famFile,as.is=T)
    
    popNames<-names(l)
    nSNP <- nrow(bim)
    nPop <- ncol(popNames)
    
    nInd <- sapply(names(l),function(x) sum(x==fam$V1))
    
    chr<-bim[,1]
    rs<-bim[,2]
    pos<-bim[,4]
    l2<-list()
    for(name in names(l)){        
        l2[[paste0(name,"MAF")]]<-l[[name]]$MAF        
    }
    freq<-do.call(cbind,l2)
    colnames(freq)<-popNames
    
} else if(filterFile!=""){
    
    qopt<-read.table(Qfile,as.is=T)
    fopt<-read.table(Ffile,as.is=T)
    beagle<-read.table(beagleFile,as.is=T,h=T)
    snpInfo<-read.table(filterFile,as.is=T,h=T)
    popInfo<-read.table(popsFile,as.is=T,h=F)
    
    ## normally do like this 
    ##popNames<-unique(popInfo$V2)
    ## ISSUE is that which column in .qopt represent which pop is kind of random
    ## min threshold or lowest max value of a row min(apply(qopt,1,max))
    ## perhaps do like this sapply(1:9,function(x) unique(popInfo[which(qopt[,x]>0.8),2]))
      
    nSNP <- nrow(snpInfo)
    nPop <- length(popNames)
    
    colnames(qopt)<-popNames
    colnames(fopt)<-popNames
    nInd <- sapply(popNames,function(x) sum(qopt[,x]))
  
    chr<-sapply(snpInfo$marker,function(x) as.numeric(gsub("chr","",unlist(strsplit(x,"_"))[1])))
    rs<-snpInfo$marker
    pos<-sapply(snpInfo$marker,function(x) as.numeric(unlist(strsplit(x,"_"))[2]))
    
 
    l2<-list()
    for(name in popNames){      
        l2[[paste0(name,"MAF")]]<-fopt[,name]        
    }
    
    freq<-do.call(cbind,l2)    
    colnames(freq)<-popNames

} else{

    qopt<-read.table(Qfile,as.is=T)
    fopt<-read.table(Ffile,as.is=T)
    bim<-read.table(bimFile,as.is=T,h=T)
    popInfo<-read.table(popsFile,as.is=T,h=F)
    
    ## normally do like this 
    ##popNames<-unique(popInfo$V2)
    ## ISSUE is that which column in .qopt represent which pop is kind of random
    ## min threshold or lowest max value of a row min(apply(qopt,1,max))
    ## perhaps do like this sapply(1:9,function(x) unique(popInfo[which(qopt[,x]>0.8),2]))
      
    nSNP <- nrow(bim)
    nPop <- length(popNames)
    
    colnames(qopt)<-popNames
    colnames(fopt)<-popNames
    nInd <- sapply(popNames,function(x) sum(qopt[,x]))
  
    chr<-bim[,1]
    rs<-bim[,2]
    pos<-bim[,4]
    
    l2<-list()
    for(name in popNames){      
        l2[[paste0(name,"MAF")]]<-fopt[,name]        
    }
    
    freq<-do.call(cbind,l2)    
    colnames(freq)<-popNames
    
}
    
    
################3


## might not need this - for if password is needed
##if(passWord!=""){
##rowNumber<-unlist(strsplit(system("grep -n '## put password here:' server.R",intern=T),":"))[1]
##system(paste0("sed '",rowNumber,"s/^/passWord<-c('",passWord,"')/' --in-place server.R"))
##}

genes<-read.table(fileGenes,h=T,as.is=T,comment.char="?")    
    
save(chr,rs,pos,file="data/snpInfo.Rdata")

save(nPop,nSNP,popNames,nInd,file="data/info.Rdata")

con<-file("data/freq.bin","wb")
gg<-as.vector(as.matrix(freq))
writeBin(gg,con)
close(con)

SNPsInChr<-table(chr)
write.table(as.matrix(SNPsInChr),"data/SNPinChr.txt",col=F,row=F,quote=F)

write.table(genes,"data/genes.anno",col=T,row=F,qu=F)

genes$chr<-sapply(genes$chrom, function(x) gsub(pattern = "chr","",x))

write.table(genes[,c("chr", "cdsStart", "cdsEnd", "name2")],"data/annoSelectedCols.txt",col=T,row=F,qu=F)

system("Rscript scripts/forGenePBS.R",intern=T)




