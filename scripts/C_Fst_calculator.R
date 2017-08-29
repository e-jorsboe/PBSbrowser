library(inline)

### ### C like R function for calculating Fst windows

calculate_Fst_for<-function(pos,vBetween,vTotal,nSNP,Vcombined,Vmidpos,regionSize,winSize){  
  
  sumVbetween<-integer(1)
  sumVtotal<-integer(1)
  minPosI=1;
  maxPosI=1
  for(i in 1:nSNP){
    
    if(pos[i]-pos[minPosI]>winSize){
      Vcombined[i]=sumVbetween/sumVtotal
      Vmidpos[i]=(pos[maxPosI]-pos[minPosI])/2+pos[minPosI]
      regionSize[i]=maxPosI-minPosI+1
      cat("region",pos[minPosI],pos[maxPosI],"Vcom",Vcombined[i],Vmidpos[i],regionSize[i],"\n")
      count=1
      #remove left side of region until pos[i] within winSize
      #pos[i] is the position of the next SNP
      while(pos[i]-pos[minPosI]>winSize){
        sumVbetween=sumVbetween-vBetween[minPosI]
        sumVtotal=sumVtotal-vTotal[minPosI]
        minPosI=minPosI+1         
      }
    }
    #expand right side
    maxPosI = i
    sumVbetween=sumVbetween+vBetween[i]
    sumVtotal=sumVtotal+vTotal[i]  
    
  }
  
  return(list(Vcombined,Vmidpos,regionSize))
}

# it calculates the sums of Vbetweem and Vtotal in each window, and then Fst=sum(Vbetween)/sum(Vtotal)
# an example of running this script is nohup Rscript C_Fst_calculator.R emilDataPH.Rdata &
# it writes a file with the Fst values (called Vcombined) 

###############################################




likeCPP_input<-signature(pos="integer",vBetween="numeric",vTotal="numeric",nSNP="integer",Vcombined="numeric",Vmidpos="numeric",regionSize="integer",winSize="numeric")

likeCPP_code<-"
double sumVbetween=0;
double sumVtotal=0;
int minPosI=0;
int maxPosI=0;
  for(int i=0;i<nSNP[0];i++){

  if(pos[i]-pos[minPosI]>winSize[0]){
    Vcombined[i]=sumVbetween/sumVtotal;
    Vmidpos[i]=(pos[maxPosI]-pos[minPosI])*1.0/2+pos[minPosI];
    regionSize[i]=maxPosI-minPosI+1;
    while(pos[i]-pos[minPosI]>winSize[0]){
      sumVbetween=sumVbetween-vBetween[minPosI];
      sumVtotal=sumVtotal-vTotal[minPosI];
      minPosI++;         
    }
  }

  maxPosI = i;
  sumVbetween=sumVbetween+vBetween[i];
  sumVtotal=sumVtotal+vTotal[i];

  }
"

fns <- cfunction(list(winSjov=likeCPP_input),
                 list(likeCPP_code),
                 convention=".C", cxxargs="-O3", cppargs="-O3",language="C++")
calculate_Fst_forCpp<-fns[["winSjov"]]


###############################

### R version too slow anyway

### essentially a windows Fst calculator, that keeps track of genomewide position
### so it is passed a variable with the number of SNPs in each chr
### enabling it to know how far in the genome it is and the it loops through each chromosome
### REMEMBER THE \\ WHEN \


likeCPP_input2<-signature(pos="integer",chr="integer",vBetween="numeric",vTotal="numeric",nSNP="integer",Vcombined="numeric",Vmidpos="numeric",regionSize="integer",winSize="numeric",SNPsinChr="integer")

likeCPP_code2<-'
int how_far=0;
/* char *str="The force is strong within me"; 
fprintf(stdout,"%s\\n",str); */
for(int j=0;j<22;j++){

double sumVbetween=0;
double sumVtotal=0;
int minPosI=0+how_far;
int maxPosI=0+how_far;
  for(int i=0;i<SNPsinChr[j];i++){

  if(pos[i+how_far]-pos[minPosI]>winSize[0]){
    Vcombined[i+how_far]=sumVbetween/sumVtotal;
    Vmidpos[i+how_far]=(pos[maxPosI]-pos[minPosI])*1.0/2+pos[minPosI];
    regionSize[i+how_far]=maxPosI-minPosI+1;
    while(pos[i+how_far]-pos[minPosI]>winSize[0]){
      sumVbetween=sumVbetween-vBetween[minPosI];
      sumVtotal=sumVtotal-vTotal[minPosI];
      minPosI++;         
    }
  }

  maxPosI = i+how_far;
  sumVbetween=sumVbetween+vBetween[i+how_far];
  sumVtotal=sumVtotal+vTotal[i+how_far];

  }
how_far+=SNPsinChr[j];

}
'

fns <- cfunction(list(winSjov=likeCPP_input2),
                 list(likeCPP_code2),
                 convention=".C", cxxargs="-O3", cppargs="-O3",language="C++")
calculate_genomeFst_forCpp<-fns[["winSjov"]]


############# Fst gene calculator



geneFstCPPinput<-signature(vBetween="numeric",vTotal="numeric",nSNP="integer",Vcombined="numeric",SNPsinGene="integer",genes="integer",backtrace="integer")


geneFstCPPcode<-'
int how_far=0;


for(int j=0;j<genes[0];j++){

double sumVbetween=0;
double sumVtotal=0;

  for(int i=0;i<SNPsinGene[j];i++){

      sumVbetween = vBetween[i+how_far] + sumVbetween;
      sumVtotal = vTotal[i+how_far] + sumVtotal;

    }
Vcombined[j] = sumVbetween / sumVtotal;
// idea is that the how_far tracking where in genome is
// put back because snp belong to more genes
how_far = SNPsinGene[j] + how_far - backtrace[j];

}
'

fns <- cfunction(list(winSjov=geneFstCPPinput),
                 list(geneFstCPPcode),
                 convention=".C", cxxargs="-O3", cppargs="-O3",language="C++")
geneFstCPP<-fns[["winSjov"]]
