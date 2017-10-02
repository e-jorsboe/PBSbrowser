library(inline)

## this is the one that should be changed!!

read1 <- function(pop,popNames1,theFile,nSNP){
  
  w <- which(popNames1==pop)
  con<-file(theFile,"rb")
  seek(con,where=nSNP*(w-1)*8)
  r<-readBin(con,"double",nSNP,size=8)
  close(con)
  r
}



## nInd,nPop,nSNP
load("data/info.Rdata")
## chr,pos
load("data/snpInfo.Rdata")

f<-list()

for(i in 1:length(popNames)){
    f[[i]]<-read1(popNames[i],popNames,"data/freq.bin",nSNP)
}

fall<-do.call(cbind,f)
colnames(fall)<-popNames

## needs calculating for table in PBS browser, closest gene and genes in window functionality
## use function closestGene in calcClosestGeneSingleSNP.R


genes<-read.table("data/genes.anno",as.is=T,head=T,comment.char="")

maxChr<-max(chr)
genes<-genes[ genes$chrom%in%sapply(1:maxChr,function(x) paste("chr",x,sep="")),]
genes$chr<-sapply(genes$chrom, function(x) unlist(strsplit(x,"chr"))[2])
genes2<-lapply(unique(genes$chr), function(x) genes[ genes$chr==x,])
names(genes2)<-unique(genes$chr)


closestGene<-function(i,genes){
  pos1<-pos[i]
  chr1<-chr[i]
  
  
  minDist<-sort.list(c(abs(pos1-genes[[chr1]][ ,"txStart"]),abs(pos1-genes[[chr1]][ ,"txEnd"])))[1] # find which distance is smallest
  minDist2<-ifelse(minDist>nrow(genes[[chr1]]),minDist-nrow(genes[[chr1]]),minDist) # because if closer to end, will be in second half of list indexes + nrow(genes)
  return(genes[[chr1]][minDist2 ,"name2"])
}

system.time(all <- parallel::mclapply(1:length(pos),function(x)closestGene(x,genes=genes2),mc.cores=10))
closestGeneSingleMarker<-cbind(chr,pos,unlist(all))
colnames(closestGeneSingleMarker)[3]<-"gene"
save(closestGeneSingleMarker,file="data/closestGeneSingleMarker.Rdata")


genesNotDup<-genes[ !duplicated(genes$name2),]
n<-length(pos)
geneMember<-logical(n)
geneNames<-character(n)


#########################################################################################

inGene<-function(i){
  
  geneMember<-genesNotDup$txStart<=pos[i] & genesNotDup$txEnd>=pos[i] & chr[i]==genesNotDup$chr
  if(any(geneMember)){
      
      geneNames<-paste(genesNotDup[geneMember,"name2"],collapse=",")
  } else{
    geneNames<-""
  }
  return(geneNames)
}

##########################################################################################################

system.time(all <- parallel::mclapply(1:length(pos), inGene,mc.cores=10))

## this has chr, pos and assigned genes plus freqs for each snp, has pos rows
## if snp in more genes then gene column is "gene1,gene2,gene3"
near<-cbind(chr,pos,unlist(all),fall)
colnames(near)<-c("chr","pos","gene",colnames(fall))
            
near[ near[,"gene"]=="","gene"]<-"NONE"

## this gives vector of one big string like this pos&chr&genes&freqs
## more strings per snp
k<-unlist(apply(near,1, function(x) sapply(unlist(strsplit(x[3],",")), function(y) paste(x[1],"&",x[2],"&",y,"&",paste(x[4:length(x)],collapse = "&"),sep=""))))
## gives list of each entry for each snp 
k2<-sapply(k, function(x) strsplit(x,"&"))
k3<-do.call(rbind,k2)
colnames(k3)<-c("chr","pos","gene",popNames)
k4<-k3
rownames(k4)<-1:nrow(k4)

assignedGenes2<-as.data.frame(k4[ k4[,"gene"]!="NONE",],stringsAsFactors=F)
uniqueAssignedGenes2<-assignedGenes2[ !duplicated(assignedGenes2[,"gene"]),]
uniqueAssignedGenes2[,"rightOrder"]<-1:nrow(uniqueAssignedGenes2)
uniqueAssignedGenes2<-uniqueAssignedGenes2[ order(uniqueAssignedGenes2$gene),]
geneCount<-table(assignedGenes2$gene)

## how many snps is in each gene, has to be ordered genomically
y<-as.data.frame(geneCount[order(uniqueAssignedGenes2$rightOrder)])
y[,2]<-rownames(y)

## for assigned genes I have to know which gene the SNP belong to
assignedGenes2$whichGeneIndex<-unlist(parallel::mclapply(1:nrow(assignedGenes2), function(x) which(assignedGenes2$gene[x]==y[,2]),mc.cores=10))

## for getting version without bloody levels

SNPsInGene<-assignedGenes2
numberGenes<-length(unique(SNPsInGene))

geneRange<-apply(as.data.frame(unique(SNPsInGene$gene)),1, function(x) genesNotDup[ genesNotDup$name2==x,c("chr","txStart","txEnd","name2")])
geneRange2<-do.call(rbind,geneRange)


## count how many SNPs in each gene
manyGene<-table(SNPsInGene$gene)
geneOrder<-cbind(as.character(unique(SNPsInGene$gene)),1:length(unique(SNPsInGene$gene)))
geneOrder2<-geneOrder[ order(geneOrder[,1]),]
manyGene2<-cbind((manyGene),geneOrder2[,2])
manyGene3<-manyGene2[ order(as.numeric(manyGene2[,2])),]
geneRange<-cbind(chr=as.numeric(geneRange2[,1]),txStart=as.numeric(geneRange2[,"txStart"]),txEnd=as.numeric(geneRange2[,"txEnd"]),
                 name2=geneRange2[,"name2"],range=paste(geneRange2[,2],"-",geneRange2[,3],sep=""),numberSNPs=as.numeric(manyGene3[,1]))


save(SNPsInGene,geneRange,file="data/geneInfo.Rdata")
