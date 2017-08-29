
source("/home/emil/native_americans_2015/C_Fst_calculator.R")

###

#load("/home/cetp/native/analyses/fst/genotypesPelChbCeu.Rdata")
#load("/home/cetp/native/analyses/fst/bimPelChbCeuFile.Rdata")
#bim<-bim[intersect(rownames(bim),rownames(geno)),]


## read in pop data, and designate pops, 1=PEL, 2=CEU, 3=CHB/HAN
#fam<-read.table("/home/cetp/native/analyses/admixture/5admPel.ind.txt")
#pop<-ifelse(fam$V2=="PEL",1,ifelse(fam$V2=="CEU",2,ifelse(fam$V2=="CHB",3,fam$V2)))


## only run on subset of data, e.g. Chromosome 11
#chr11<-geno[bim$V1=="11",]
#genot<-t(chr11)
#bimt<-bim[bim$V1=="11",]

## pel and ceu 1,2
## pel and han 1,3 
## han and ceu 2,3 
## does it all in one go

# this calculates the Nei Fst values, or the upper and lower part of the Nei eq
# which it writes to a file called nei_Fst_$pop1,$pop2


#######

reynolds_matrix<-function(geno,populations,bim,calc_pbs){
  ### remove missing data from it 
  
  ### so first calculate Freq for each population
  
  # removing rows with a NA
  geno<-geno[,apply(geno,2,function(x) all(!is.na(x)))]
  
  # for the out SNP ids
  ids<-colnames(geno)
  pos<-sapply(ids,function(x) unlist(strsplit(x,"_"))[2])
  
  bim<-bim[bim$V4%in%pos,]
  
  rownames(geno)<-populations
  
  # freqs counting the number of minor alleles and divinding by 2*n 
  f<-sapply(sort(unique(populations)), function(x) apply(geno[rownames(geno) %in% x,],2,function(y) sum(y) / as.numeric(length(y)*2)))
  
  # if want to try and emulate Matteo's student
  #f["11_61597212",2]<-0.05 CEU of 5 %
  #f["11_61597212",3]<-0.25 HAN of 25 %
  
  q<-1-f
  alpha <- 1 - (f**2 + q**2)
  
  # getting the combinations for the populations to calculate Fst between
  combinations<-unlist(sapply(1:(max(populations)-1),function(x) paste(x,unique(populations[populations>x]),sep = ",")))

  # do function instead and put in sapply
  
  calc_al<-function(x){
    p1<-f[,as.numeric(unlist(strsplit(x,","))[1])]
    p2<-f[,as.numeric(unlist(strsplit(x,","))[2])]**2  # (p1-p2)**2
    q1<-q[,as.numeric(unlist(strsplit(x,","))[1])]
    q2<-q[,as.numeric(unlist(strsplit(x,","))[2])]**2 # (q1-q2)**2
    N1<-length(populations[populations==as.numeric(unlist(strsplit(x,","))[1])])
    N2<-length(populations[populations==as.numeric(unlist(strsplit(x,","))[2])])  # N1+N2
    alpha1<-alpha[,as.numeric(unlist(strsplit(x,","))[1])]
    alpha2<-alpha[,as.numeric(unlist(strsplit(x,","))[2])]
    #print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    # do formula with those variables like 
    al2 <- 1/2 * ( (p1-p2)^2 + (q1-q2)^2) - (N1+N2) * (N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
    return(al2)
    #return((0.5*((p1-p2)**2 + (q1-q2)**2) - (N1+N2) * (N1*alpha1 + N2*alpha2)) / (4*N1*N2*(N1+N2-1)))
  }
  
  
  al<-sapply(combinations, calc_al)

 
  calc_bal<-function(x){
    p1<-f[,as.numeric(unlist(strsplit(x,","))[1])]
    p2<-f[,as.numeric(unlist(strsplit(x,","))[2])]**2  # (p1-p2)**2
    q1<-q[,as.numeric(unlist(strsplit(x,","))[1])]
    q2<-q[,as.numeric(unlist(strsplit(x,","))[2])]**2  # (q1-q2)**2
    N1<-length(populations[populations==as.numeric(unlist(strsplit(x,","))[1])])
    N2<-length(populations[populations==as.numeric(unlist(strsplit(x,","))[2])])  # N1+N2
    alpha1<-alpha[,as.numeric(unlist(strsplit(x,","))[1])]
    alpha2<-alpha[,as.numeric(unlist(strsplit(x,","))[2])]
    # do formula with those variables like 
    
    #return((0.5*((p1-p2)**2 + (q1-q2)**2) + (4*N1*N2-N1-N2)*(N1*alpha1 + N2*alpha2)) / (4*N1*N2*(N1+N2-1)))
    bal2 <- 1/2 * ( (p1-p2)^2 + (q1-q2)^2) + (4*N1*N2-N1-N2)*(N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
    return(bal2)
    
  }
  
  bal<-sapply(combinations,calc_bal)

  if(calc_pbs==0) {
    # return Fst values top and bottom
    write.table(cbind(bim$V1,bim$V4,al[,1],bal[,1]),"unadmixed_Rey_Fst_1,2",col.names = c("chr","pos","rey_top","rey_bottom"),row=F,quote=F)
    write.table(cbind(bim$V1,bim$V4,al[,2],bal[,2]),"unadmixed_Rey_Fst_1,3",col.names = c("chr","pos","rey_top","rey_bottom"),row=F,quote=F)
    write.table(cbind(bim$V1,bim$V4,al[,3],bal[,3]),"unadmixed_Rey_Fst_2,3",col.names = c("chr","pos","rey_top","rey_bottom"),row=F,quote=F)
  } else{
    # for a single locus
  
  # for a single locus
    fst<- (al / bal)
    
    pbs<-apply(fst,1,function(x) (-log(1-x[1]) + -log(1-x[2]) - -log(1-x[3])) / 2)
    output<-cbind(bim$V1,bim$V4,pbs)
    colnames(output)<-c("chr","pos","PBS")
    write.table(output,paste("chr11_PBS_PELCEUCHB_unadmixfreq.txt",sep=""),col=T,row=F,quote=F)
    return(pbs)
  }
}

# give vector of pop size of each population

# give freqs directly ### HOW TO KNOW IF DAF OR AAF?


#admix_indis<-read.table("/pontus/data/cetp/native/analyses/admixture/PEL/merged.3.Q_17",as.is=T)
#colnames(admix_indis)<-c("CHB","CEU", "PEL")
#admix_loci<-read.table("/pontus/data/cetp/native/analyses/admixture/PEL/merged.3.P_17.gz",as.is=T)

#fam<-read.table("/pontus/data/cetp/native/raw/plink/CEUCHBPEL/merged.fam",as.is=T)
#bim<-read.table("/pontus/data/cetp/native/raw/plink/CEUCHBPEL/merged.bim",as.is=T)

#chr11_admix_loci<-admix_loci[which(bim$V1==11),]

#pop<-ifelse(fam$V2=="PEL",1,ifelse(fam$V2=="CEU",2,ifelse(fam$V2=="CHB",3,fam$V2)))
#pop1<-as.numeric(pop)

reynolds_admix<-function(freq,populations,pop_sizes,bim,orderPBS,windows){
  
  #freq=fall[right_pos,];populations=c(1,2,3);pop_sizes=nInd;bim=chr_pos[right_pos,];orderPBS = c("PEL","CEU","CHB");windows=0;
  # supply pop1, pop2 and outgroup as arguments
  # and then [pop1,pop2], [pop1,pop3] & [pop2,pop3] in fst matrix for when doing the PBS
  
  #freq<-cbind(freq[,orderPBS[1]],freq[,orderPBS[2]],freq[,orderPBS[3]])
  #colnames(freq)<-c(orderPBS[1],orderPBS[2],orderPBS[3])
  #freq = chr11_admix_loci
  #populations = pop1
  #pop_sizes = c(sum(admix_indis$PEL),sum(admix_indis$CEU),sum(admix_indis$CHB))
  #bim = bim[bim$V1==11,]
  # order populations differently so that 1=PEL,2=CEU,3=CHB
  #freq<-fall
  #pop_sizes<-nInd
  # if want to try and emulate Matteo's student
  #f["11_61597212",2]<-0.05 CEU of 5 %
  #f["11_61597212",3]<-0.25 HAN of 25 %
  
  # rearrange admixture output so 1=PEL,2=CEU,3=CHB
  #freq1<-cbind(freq[,3],freq[,2],freq[,1])
  #pop_sizes1<-c(pop_sizes[3],pop_sizes[2],pop_sizes[1])
  
  # print(head(freq))
  q<-1-freq
  alpha <- 1 - (freq**2 + q**2)
  
  # getting the combinations for the populations to calculate Fst between
  #combinations2<-unlist(sapply(1:(max(populations)-1),function(x) paste(x,unique(populations[populations>x]),sep = ",")))
  
  #combinations<-unlist(sapply(orderPBS,function(x) paste(x,popNames[popNames!=x],sep = ",")))
  
  combinations<-c()
  combinations[1]<-paste(orderPBS[1],orderPBS[2],sep = ",")
  combinations[2]<-paste(orderPBS[1],orderPBS[3],sep = ",")
  combinations[3]<-paste(orderPBS[2],orderPBS[3],sep = ",")
  
  
  calc_al<-function(f1,f2,N1,N2){
    q1<-1-f1
    q2<-1-f2
    alpha1<-1 - (f1**2 + q1**2)
    alpha2<-1 - (f2**2 + q2**2)
    #print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    # do formula with those variables like 
    al2 <- 1/2 * ( (f1-f2)^2 + (q1-q2)^2) - (N1+N2) * (N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
    return(al2)
    #return((0.5*((p1-p2)**2 + (q1-q2)**2) - (N1+N2) * (N1*alpha1 + N2*alpha2)) / (4*N1*N2*(N1+N2-1)))
  }
  
  #calc_al(fall[right_pos,"CEU"],fall[right_pos,"PEL"],nInd["CEU"],nInd["PEL"])/calc_bal(fall[right_pos,"CEU"],fall[right_pos,"PEL"],nInd["CEU"],nInd["PEL"])
  
  al<-cbind(calc_al(freq[,orderPBS[1]],freq[,orderPBS[2]],pop_sizes[orderPBS[1]],pop_sizes[orderPBS[2]]),
             calc_al(freq[,orderPBS[1]],freq[,orderPBS[3]],pop_sizes[orderPBS[1]],pop_sizes[orderPBS[3]]),
             calc_al(freq[,orderPBS[2]],freq[,orderPBS[3]],pop_sizes[orderPBS[2]],pop_sizes[orderPBS[3]]))
  
  calc_al_old<-function(x){
    
    p1<-freq[,unlist(strsplit(x,","))[1]]
    p2<-freq[,unlist(strsplit(x,","))[2]]  # (p1-p2)**2
    q1<-q[,unlist(strsplit(x,","))[1]]
    q2<-q[,unlist(strsplit(x,","))[2]] # (q1-q2)**2
    N1<-pop_sizes[unlist(strsplit(x,","))[1]]
    N2<-pop_sizes[unlist(strsplit(x,","))[2]]
    alpha1<-alpha[,unlist(strsplit(x,","))[1]]
    alpha2<-alpha[,unlist(strsplit(x,","))[2]]
    #print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    # do formula with those variables like 
    al2 <- 1/2 * ( (p1-p2)^2 + (q1-q2)^2) - (N1+N2) * (N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
    return(al2)
    #return((0.5*((p1-p2)**2 + (q1-q2)**2) - (N1+N2) * (N1*alpha1 + N2*alpha2)) / (4*N1*N2*(N1+N2-1)))
  }
  
  calc_al2<-function(x){
    
    p1<-freq[,as.numeric(unlist(strsplit(x,","))[1])]
    p2<-freq[,as.numeric(unlist(strsplit(x,","))[2])]  # (p1-p2)**2
    q1<-q[,as.numeric(unlist(strsplit(x,","))[1])]
    q2<-q[,as.numeric(unlist(strsplit(x,","))[2])] # (q1-q2)**2
    N1<-pop_sizes[as.numeric(unlist(strsplit(x,","))[1])]
    N2<-pop_sizes[as.numeric(unlist(strsplit(x,","))[2])]
    alpha1<-alpha[,as.numeric(unlist(strsplit(x,","))[1])]
    alpha2<-alpha[,as.numeric(unlist(strsplit(x,","))[2])]
    #print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    # do formula with those variables like 
    al2 <- 1/2 * ( (p1-p2)^2 + (q1-q2)^2) - (N1+N2) * (N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
    return(al2)
    #return((0.5*((p1-p2)**2 + (q1-q2)**2) - (N1+N2) * (N1*alpha1 + N2*alpha2)) / (4*N1*N2*(N1+N2-1)))
  }
  
  #al<-sapply(combinations, calc_al)
  #al2<-sapply(combinations2, calc_al2)
  #colnames(al)<-unlist(combinations)
  
  calc_bal<-function(f1,f2,N1,N2){
    q1<-1-f1
    q2<-1-f2
    alpha1<-1 - (f1**2 + q1**2)
    alpha2<-1 - (f2**2 + q2**2)
    #print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    # do formula with those variables like 
    bal2 <- 1/2 * ( (f1-f2)^2 + (q1-q2)^2) + (4*N1*N2-N1-N2)*(N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
    return(bal2)
    #return((0.5*((p1-p2)**2 + (q1-q2)**2) - (N1+N2) * (N1*alpha1 + N2*alpha2)) / (4*N1*N2*(N1+N2-1)))
  }
  
  calc_bal_old<-function(x){
    p1<-freq[,unlist(strsplit(x,","))[1]]
    p2<-freq[,unlist(strsplit(x,","))[2]]  # (p1-p2)**2
    q1<-q[,unlist(strsplit(x,","))[1]]
    q2<-q[,unlist(strsplit(x,","))[2]] # (q1-q2)**2
    N1<-pop_sizes[unlist(strsplit(x,","))[1]]
    N2<-pop_sizes[unlist(strsplit(x,","))[2]]
    alpha1<-alpha[,unlist(strsplit(x,","))[1]]
    alpha2<-alpha[,unlist(strsplit(x,","))[2]]
    #print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    # do formula with those variables like 
    bal2 <- 1/2 * ( (p1-p2)^2 + (q1-q2)^2) + (4*N1*N2-N1-N2)*(N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
    return(bal2)
    #return((0.5*((p1-p2)**2 + (q1-q2)**2) - (N1+N2) * (N1*alpha1 + N2*alpha2)) / (4*N1*N2*(N1+N2-1)))
  }
  
  
  calc_bal2<-function(x){
    print(as.numeric(unlist(strsplit(x,","))[1]))
    p1<-freq[,as.numeric(unlist(strsplit(x,","))[1])]
    p2<-freq[,as.numeric(unlist(strsplit(x,","))[2])]  # (p1-p2)**2
    q1<-q[,as.numeric(unlist(strsplit(x,","))[1])]
    q2<-q[,as.numeric(unlist(strsplit(x,","))[2])] # (q1-q2)**2
    N1<-pop_sizes[as.numeric(unlist(strsplit(x,","))[1])]
    N2<-pop_sizes[as.numeric(unlist(strsplit(x,","))[2])]
    alpha1<-alpha[,as.numeric(unlist(strsplit(x,","))[1])]
    alpha2<-alpha[,as.numeric(unlist(strsplit(x,","))[2])]
    #print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    # do formula with those variables like 
    bal2 <- 1/2 * ( (p1-p2)^2 + (q1-q2)^2) + (4*N1*N2-N1-N2)*(N1*alpha1 + N2*alpha2) / (4*N1*N2*(N1+N2-1))
    return(bal2)
    #return((0.5*((p1-p2)**2 + (q1-q2)**2) - (N1+N2) * (N1*alpha1 + N2*alpha2)) / (4*N1*N2*(N1+N2-1)))
  }
  
  bal<-cbind(calc_bal(freq[,orderPBS[1]],freq[,orderPBS[2]],pop_sizes[orderPBS[1]],pop_sizes[orderPBS[2]]),
             calc_bal(freq[,orderPBS[1]],freq[,orderPBS[3]],pop_sizes[orderPBS[1]],pop_sizes[orderPBS[3]]),
             calc_bal(freq[,orderPBS[2]],freq[,orderPBS[3]],pop_sizes[orderPBS[2]],pop_sizes[orderPBS[3]]))
  
  # looks like the values have just been changed around but otherwise the same, MIGHT BE THE PROBLEM
  #bal<-sapply(combinations,calc_bal)
  #bal2<-sapply(combinations2,calc_bal2)
  #colnames(bal)<-unlist(combinations)
  
  #print(tail(head(al/bal,n=3050),n=1))
  
  #index12<-which(grepl(orderPBS[1],combinations) &  grepl(orderPBS[2],combinations))
  #index13<-which(grepl(orderPBS[1],combinations) &  grepl(orderPBS[3],combinations))
  #index23<-which(grepl(orderPBS[2],combinations) &  grepl(orderPBS[3],combinations))
  
  if(windows==0){
  # for a single locus
    fst<- (al / bal)
 
 
    #pbs<-apply((al/bal),1,function(x) (-log(1-x[index12]) + -log(1-x[index13]) - -log(1-x[index23])) / 2)
    #pbs2<-apply((al2/bal2),1,function(x) (-log(1-x[index12]) + -log(1-x[index13]) - -log(1-x[index23])) / 2)
    pbs<-apply(fst,1,function(x) (-log(1-x[1]) + -log(1-x[2]) - -log(1-x[3])) / 2)
    output<-cbind(bim$V1,bim$V4,pbs,fst[,1],fst[,2],fst[,3])
    colnames(output)<-c("chr","pos","PBS",paste("Fst",orderPBS[1],orderPBS[2],sep=""),paste("Fst",orderPBS[1],orderPBS[3],sep=""),paste("Fst",orderPBS[2],orderPBS[3],sep=""))
    #write.table(output,paste("chr11_PBS_PELCEUCHB_admixfreq.txt",sep=""),col=T,row=F,quote=F)
    return(output)
  } else {   
      winFst<-list()
      # calc Fst for windows with C_Fst_calculator
      for(p in 1:ncol(al)){
        
          nSNP<-nrow(bim)
          Vcombined<-numeric(nSNP)
          regionSize<-integer(nSNP)
          Vmidpos<-numeric(nSNP)
          pos=as.numeric(bim$V4)
          vBetween=as.numeric(al[,p])
          vTotal=as.numeric(bal[,p])
          winSize=as.numeric(windows)
          maxChr<-as.numeric(max(bim$V1))
          
          cal2<-calculate_Fst_forCpp(pos,vBetween,vTotal,nSNP=nSNP,Vcombined=Vcombined,Vmidpos=Vmidpos,winSize=winSize,regionSize=regionSize,maxChr=maxChr)
          res<-do.call(cbind,cal2)
          res2<-res[res[,7]>0,]
          winFst[[p]]<-res2[,"Vcombined"]
          names(winFst)[p]<-colnames(al)[p]
          win_midpos<-res2[,"Vmidpos"]   
          region_size<-res2[,"regionSize"]
      }   
      
      #pbs<-sapply(1:length(winFst[[1]]),function(x) (-log(1-winFst[[index12]][x]) + -log(1-winFst[[index13]][x]) - -log(1-winFst[[index23]][x])) / 2)
      pbs<-sapply(1:length(winFst[[1]]),function(x) (-log(1-winFst[[1]][x]) + -log(1-winFst[[2]][x]) - -log(1-winFst[[3]][x])) / 2)
      pbs2<-cbind(unique(bim$V1),win_midpos,pbs,region_size,winFst[[1]],winFst[[2]],winFst[[3]])
      colnames(pbs2)<-c("chr","pos","PBS","SNPsinWin",paste("Fst",orderPBS[1],orderPBS[2],sep=""),paste("Fst",orderPBS[1],orderPBS[3],sep=""),paste("Fst",orderPBS[2],orderPBS[3],sep=""))
      #write.table(cbind(bim$V1,bim$V4,al[,1],bal[,1]),"admixed_Rey_Fst_1,2",col.names = c("chr","pos","rey_top","rey_bottom"),row=F,quote=F)
      #write.table(cbind(bim$V1,bim$V4,al[,2],bal[,2]),"admixed_Rey_Fst_1,3",col.names = c("chr","pos","rey_top","rey_bottom"),row=F,quote=F)
      #write.table(cbind(bim$V1,bim$V4,al[,3],bal[,3]),"admixed_Rey_Fst_2,3",col.names = c("chr","pos","rey_top","rey_bottom"),row=F,quote=F)
      return(pbs2)
    }
}

#right_pos<-which(pos>=60e6 & pos<=62e6 & chr==11)
#pbs<-reynolds_admix(freq=fall[right_pos,],populations=c(1,2,3),pop_sizes=nInd,bim=chr_pos[right_pos,],orderPBS = c("CEU","PEL","CHB"),windows=0)

#right_pos<-which(pos>=60e6 & pos<=62e6 & chr==11)
#pbs<-reynolds_admix(freq=fall[right_pos,],populations=c(1,2,3),pop_sizes=nInd,bim=chr_pos[right_pos,],orderPBS = c("CEU","PEL","CHB"),windows=0)
#pbs2<-reynolds_admix(freq=fall[right_pos,],populations=c(1,2,3),pop_sizes=nInd,bim=chr_pos[right_pos,],orderPBS = c("PEL","CEU","CHB"),windows=1)


