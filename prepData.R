arguments<-commandArgs(trailingOnly=T)

if(length(arguments)==0){
    print("Please supply arguments:")
    print("1. Folder of .frq files, folder ONLY with .frq files that should be used, names of .frq files should be same as family ID in .fam file.")
    print("2. bim file of data .frq files are based on.")
    print("3. fam file of data .frq files are based on, family ID (1st column), should be same as name of .frq files")
    print("4. gene annotation file, like ???")
    q()
          

    }

## Have to be in folder when running!
setwd(".")


system("mkdir -p data",intern=T)


## reads all .frq files in this directory
freqDirectory<-arguments[1]
bimFile<-arguments[2]
famFile<-arguments[3]
geneFile<-arguments[4]
folderName<-arguments[5]

####################################3
## files should have same name as Family ID in .fam file!!!!!!!!!!!!
########################################

## freqDirectory<-"/home/emil/waterbucks/data/generated/byChr/duplicatesMergedNind50Maf005/calledGenos/R209andGprobs095/byCountry/forShiny/" 
## bimFile<-"/home/emil/waterbucks/data/generated/byChr/duplicatesMergedNind50Maf005/calledGenos/R209andGprobs095/gooodMapabilityChrALLR209gprobs095updatedIDV3_noKOB_Geno005.bim"
## famFile<-"/home/emil/waterbucks/data/generated/byChr/duplicatesMergedNind50Maf005/calledGenos/R209andGprobs095/gooodMapabilityChrALLR209gprobs095updatedIDV3_noKOB_Geno005.fam"
## geneFile<-"/home/emil/shiny/waterbuckPBSV2/Bostau8_UMD3.1.1_anno.gz"
## folderName<-"genericPBScalculator"

write(folderName,"folderName")

freqFiles<-list.files(freqDirectory,full=T,pattern=".frq$")

l<-list()
for(file in freqFiles){
    
    name<-unlist(strsplit(tail(unlist(strsplit(file,"/")),1),".frq"))[1]
    l[[name]]<-read.table(file,h=T,as.is=T)
    print("done:")
    print(file)
    
}

bim<-read.table(bimFile,as.is=T,colC=c("integer","character")[c(1,2,1,1,2,2)])
fam<-read.table(famFile,as.is=T)

genes<-read.table(geneFile,h=T,as.is=T)

popNames<-names(l)
nSNP <- nrow(bim)
nPop <- ncol(popNames)

nInd <- sapply(names(l),function(x) sum(x==fam$V1))

################3

chr<-bim[,1]
rs<-bim[,2]
pos<-bim[,4]

save(chr,rs,pos,file="data/snpInfo.Rdata")

save(nPop,nSNP,popNames,nInd,file="data/info.Rdata")

l2<-list()
for(name in names(l)){

    l2[[paste0(name,"MAF")]]<-l[[name]]$MAF

}

freq<-do.call(cbind,l2)

colnames(freq)<-popNames

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




