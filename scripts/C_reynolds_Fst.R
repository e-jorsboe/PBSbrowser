
#source("/home/emil/native_americans_2015/C_Fst_calculator.R")

## in order to generalize we would need to supply a pops*sites long vector
## same goes for popSizes, but you would need a k with number of pops
## N is a pops*sites long vector, but with different counts if missing data
## put these statements into more like functions - stop waste of code
## keep everything in vectors and use for loops - do not use apply - then C++ will come easily

C_like_reynold<-function(f1,f2,f3,N1,N2,N3,pos,chr,windows=0){
  q1<-numeric(length(f1))
  q2<-numeric(length(f1))
  q3<-numeric(length(f1))
  alpha1<-numeric(length(f1))
  alpha2<-numeric(length(f1))
  alpha3<-numeric(length(f1))
  al12<-numeric(length(f1))
  al13<-numeric(length(f1))
  al23<-numeric(length(f1))
  bal12<-numeric(length(f1))
  bal13<-numeric(length(f1))
  bal23<-numeric(length(f1))
  pbs<-numeric(length(f1))
  
  for(i in 1:length(f1)){
    q1[i]<-1-f1[i]
    q2[i]<-1-f2[i]
    q3[i]<-1-f3[i]
    alpha1[i]<-1 - (f1[i]**2 + q1[i]**2)
    alpha2[i]<-1 - (f2[i]**2 + q2[i]**2)
    alpha3[i]<-1 - (f3[i]**2 + q3[i]**2)
    #print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    # do formula with those variables like 
    # Fst 1,2
    al12[i] <- 1/2 * ( (f1[i]-f2[i])^2 + (q1[i]-q2[i])^2) - (N1+N2) * (N1*alpha1[i] + N2*alpha2[i]) / (4*N1*N2*(N1+N2-1))
    bal12[i] <- 1/2 * ( (f1[i]-f2[i])^2 + (q1[i]-q2[i])^2) + (4*N1*N2-N1-N2)*(N1*alpha1[i] + N2*alpha2[i]) / (4*N1*N2*(N1+N2-1))
    # Fst 1,3
    al13[i] <- 1/2 * ( (f1[i]-f3[i])^2 + (q1[i]-q3[i])^2) - (N1+N3) * (N1*alpha1[i] + N3*alpha3[i]) / (4*N1*N3*(N1+N3-1))
    bal13[i] <- 1/2 * ( (f1[i]-f3[i])^2 + (q1[i]-q3[i])^2) + (4*N1*N3-N1-N3)*(N1*alpha1[i] + N3*alpha3[i]) / (4*N1*N3*(N1+N3-1))
    # Fst 2,3
    al23[i] <- 1/2 * ( (f2[i]-f3[i])^2 + (q2[i]-q3[i])^2) - (N2+N3) * (N2*alpha2[i] + N3*alpha3[i]) / (4*N2*N3*(N2+N3-1))
    bal23[i] <- 1/2 * ( (f2[i]-f3[i])^2 + (q2[i]-q3[i])^2) + (4*N2*N3-N2-N3)*(N2*alpha2[i] + N3*alpha3[i]) / (4*N2*N3*(N2+N3-1))
  }
  if(windows==1){
    return(c(al12,bal12,al13,bal13,al23,bal23))
  } else{
    for(i in (1:length(f1))){
      pbs[i]<- (-log(1-(al12[i]/bal12[i])) + -log(1-(al13[i]/bal13[i])) - -log(1-(al23[i]/bal23[i]))) / 2
      
    }
    return(cbind(chr,pos,pbs))
  }
}




C_like_nei<-function(f1,f2,f3,pos,chr,windows=0){
  al12<-numeric(length(f1))
  al13<-numeric(length(f1))
  al23<-numeric(length(f1))
  bal12<-numeric(length(f1))
  bal13<-numeric(length(f1))
  bal23<-numeric(length(f1))
  pbs<-numeric(length(f1))
  
  for(i in 1:length(f1)){
    #print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    # do formula with those variables like 
    # Fst 1,2
    al12[i] <- (f1[i]+f2[i])*(f1[i]+f2[i])
    bal12[i] <- 2*((f1[i]+f2[i])/2)*(1-((f1[i]+f2[i])/2))
    # Fst 1,3
    al13[i] <- (f1[i]+f3[i])*(f1[i]+f3[i]) 
    bal13[i] <- 2*((f1[i]+f3[i])/2)*(1-((f1[i]+f3[i])/2))
    # Fst 2,3
    al23[i] <- (f2[i]+f3[i])*(f2[i]+f3[i]) 
    bal23[i] <- 2*((f2[i]+f3[i])/2)*(1-((f2[i]+f3[i])/2))
  }
  if(windows==1){
    return(c(al12,bal12,al13,bal13,al23,bal23))
  } else{
    for(i in (1:length(f1))){
      pbs[i]<- (-log(1-(al12[i]/bal12[i])) + -log(1-(al13[i]/bal13[i])) - -log(1-(al23[i]/bal23[i]))) / 2
      
    }
    return(cbind(chr,pos,pbs))
  }
}




C_like_hudson<-function(f1,f2,f3,N1,N2,N3,pos,chr,windows=0){
  q1<-numeric(length(f1))
  q2<-numeric(length(f1))
  q3<-numeric(length(f1))

  al12<-numeric(length(f1))
  al13<-numeric(length(f1))
  al23<-numeric(length(f1))
  bal12<-numeric(length(f1))
  bal13<-numeric(length(f1))
  bal23<-numeric(length(f1))
  pbs<-numeric(length(f1))
  
  for(i in 1:length(f1)){
    q1[i]<-1-f1[i]
    q2[i]<-1-f2[i]
    q3[i]<-1-f3[i]
    ##print(paste(c(p1,p2,q1,q2,N1,N2,alpha1,alpha2),collapse = "-"))
    ## do formula with those variables like 
    ## Fst 1,2
    al12[i] <- (f1[i]-f2[i])**2 - ((f1[i]*q1[i])/(N1-1)) - ((f2[i]*q2[i])/(N2-1))
    bal12[i] <- f1[i]*q2[i] + f2[i]*q1[i] 
    ## Fst 1,3
    al13[i] <- (f1[i]-f3[i])**2 - ((f1[i]*q1[i])/(N1-1)) - ((f3[i]*q3[i])/(N3-1))
    bal13[i] <- f1[i]*q3[i] + f3[i]*q1[i] 
    ## Fst 2,3
    al23[i] <- (f2[i]-f3[i])**2 - ((f2[i]*q2[i])/(N2-1)) - ((f3[i]*q3[i])/(N3-1))
    bal23[i] <- f2[i]*q3[i] + f3[i]*q2[i]
    
  }
  if(windows==1){
    return(c(al12,bal12,al13,bal13,al23,bal23))
  } else{
    for(i in (1:length(f1))){
      pbs[i]<- (-log(1-(al12[i]/bal12[i])) + -log(1-(al13[i]/bal13[i])) - -log(1-(al23[i]/bal23[i]))) / 2
      
    }
    return(cbind(chr,pos,pbs))
  }
}






#C_like_reynolds<-function(f1,f2,f3,N1,N2,N3,pos,chr)

#right_pos<-which(pos>=60000000 & pos<=62000000 & chr==11)
#right_pos<-which(chr==11)

#noWin<-C_like_reynolds(fall[right_pos,"PEL"],fall[right_pos,"CEU"],fall[right_pos,"CHB"],nInd["PEL"],nInd["CEU"],nInd["CHB"],pos[right_pos],chr[right_pos])

#right_pos<-which(pos>=60e6 & pos<=62e6 & chr==11)
#pbs<-reynolds_admix(freq=fall[right_pos,],populations=c(1,2,3),pop_sizes=nInd,bim=chr_pos[right_pos,],orderPBS = c("CEU","PEL","CHB"),windows=0)

#right_pos<-which(pos>=60e6 & pos<=62e6 & chr==11)
#pbs<-reynolds_admix(freq=fall[right_pos,],populations=c(1,2,3),pop_sizes=nInd,bim=chr_pos[right_pos,],orderPBS = c("CEU","PEL","CHB"),windows=0)
#pbs2<-reynolds_admix(freq=fall[right_pos,],populations=c(1,2,3),pop_sizes=nInd,bim=chr_pos[right_pos,],orderPBS = c("PEL","CEU","CHB"),windows=1)


################ C version


library(inline)


## both calculates PBS and the variane between populations (al) and the total (bal) variance (between + within)
likeCPP_input_ReynoldFst<-signature(f1="numeric",f2="numeric",f3="numeric",N1="numeric",N2="numeric",N3="numeric",pos="integer",chr="integer",n="integer",al12="numeric",
                         al13="numeric",al23="numeric",bal12="numeric",bal13="numeric",bal23="numeric",pbs="numeric")
likeCPP_code_ReynoldFst<-"
double* q1 = new double[n[0]];
double* q2 = new double[n[0]];
double* q3 = new double[n[0]];

double* alpha1 = new double[n[0]];
double* alpha2 = new double[n[0]];
double* alpha3 = new double[n[0]];

/* q1 double[n[0]] = { 0.0 };
q2 double[n[0]] = { 0.0 };
q3 double[n[0]] = { 0.0 };

alpha1 double[n[0]] = { 0.0 };
alpha2 double[n[0]] = { 0.0 };
alpha3 double[n[0]] = { 0.0 }; */

for(int i=0;i<n[0];i++){
   
   q1[i]=1-f1[i];
   q2[i]=1-f2[i];
   q3[i]=1-f3[i];
   alpha1[i]=1 - (f1[i]*f1[i] + q1[i]*q1[i]);
   alpha2[i]=1 - (f2[i]*f2[i] + q2[i]*q2[i]);
   alpha3[i]=1 - (f3[i]*f3[i] + q3[i]*q3[i]);
   al12[i] = 0.5*((f1[i]-f2[i])*(f1[i]-f2[i]) + (q1[i]-q2[i])*(q1[i]-q2[i])) - (N1[0]+N2[0]) * (N1[0]*alpha1[i] + N2[0]*alpha2[i]) / (4*N1[0]*N2[0]*(N1[0]+N2[0]-1));
   bal12[i] = 0.5*((f1[i]-f2[i])*(f1[i]-f2[i]) + (q1[i]-q2[i])*(q1[i]-q2[i])) + (4*N1[0]*N2[0]-N1[0]-N2[0])*(N1[0]*alpha1[i] + N2[0]*alpha2[i]) / (4*N1[0]*N2[0]*(N1[0]+N2[0]-1));
   
   al13[i] = 0.5*((f1[i]-f3[i])*(f1[i]-f3[i]) + (q1[i]-q3[i])*(q1[i]-q3[i])) - (N1[0]+N3[0]) * (N1[0]*alpha1[i] + N3[0]*alpha3[i]) / (4*N1[0]*N3[0]*(N1[0]+N3[0]-1));
   bal13[i] = 0.5*((f1[i]-f3[i])*(f1[i]-f3[i]) + (q1[i]-q3[i])*(q1[i]-q3[i])) + (4*N1[0]*N3[0]-N1[0]-N3[0])*(N1[0]*alpha1[i] + N3[0]*alpha3[i]) / (4*N1[0]*N3[0]*(N1[0]+N3[0]-1));
   
   al23[i] = 0.5*((f2[i]-f3[i])*(f2[i]-f3[i]) + (q2[i]-q3[i])*(q2[i]-q3[i])) - (N2[0]+N3[0]) * (N2[0]*alpha2[i] + N3[0]*alpha3[i]) / (4*N2[0]*N3[0]*(N2[0]+N3[0]-1));
   bal23[i] = 0.5*((f2[i]-f3[i])*(f2[i]-f3[i]) + (q2[i]-q3[i])*(q2[i]-q3[i])) + (4*N2[0]*N3[0]-N2[0]-N3[0])*(N2[0]*alpha2[i] + N3[0]*alpha3[i]) / (4*N2[0]*N3[0]*(N2[0]+N3[0]-1));
  } 
for(int i=0;i<n[0];i++){
  pbs[i]= (-log(1-(al12[i]/bal12[i])) + -log(1-(al13[i]/bal13[i])) - -log(1-(al23[i]/bal23[i]))) / 2.0;
      
}
"

fns <- cfunction(list(pbsCalculator_ReynoldFst=likeCPP_input_ReynoldFst),
                 list(likeCPP_code_ReynoldFst),
                 convention=".C", cxxargs="-O3", cppargs="-O3",language="C++")
pbsCalculator_forCpp_ReynoldFst<-fns[["pbsCalculator_ReynoldFst"]]






## both calculates PBS and the variane between populations (al) and the total (bal) variance (between + within)
likeCPP_input_NeiFst<-signature(f1="numeric",f2="numeric",f3="numeric",pos="integer",chr="integer",n="integer",al12="numeric",al13="numeric",al23="numeric",bal12="numeric",bal13="numeric",bal23="numeric",pbs="numeric")
likeCPP_code_NeiFst<-"

for(int i=0;i<n[0];i++){

   al12[i] = (f1[i]-f2[i])*(f1[i]-f2[i]);
   bal12[i] = 2*((f1[i]+f2[i])/2)*(1-((f1[i]+f2[i])/2));
   
   al13[i] = (f1[i]-f3[i])*(f1[i]-f3[i]);
   bal13[i] = 2*((f1[i]+f3[i])/2)*(1-((f1[i]+f3[i])/2));
   
   al23[i] = (f2[i]-f3[i])*(f2[i]-f3[i]);
   bal23[i] = 2*((f2[i]+f3[i])/2)*(1-((f2[i]+f3[i])/2));
  } 
for(int i=0;i<n[0];i++){
  pbs[i]= (-log(1-(al12[i]/bal12[i])) + -log(1-(al13[i]/bal13[i])) - -log(1-(al23[i]/bal23[i]))) / 2.0;
      
}
"

fns <- cfunction(list(pbsCalculator_NeiFst=likeCPP_input_NeiFst),
                 list(likeCPP_code_NeiFst),
                 convention=".C", cxxargs="-O3", cppargs="-O3",language="C++")
pbsCalculator_forCpp_NeiFst<-fns[["pbsCalculator_NeiFst"]]




####################################

## based on Hudson's Fst formula proposed in Gaurav Bhatia et al., 2013 - Genome Research
## both calculates PBS and the variane between populations (al) and the total (bal) variance (between + within)
likeCPP_input_HudsonFst<-signature(f1="numeric",f2="numeric",f3="numeric",N1="numeric",N2="numeric",N3="numeric",pos="integer",chr="integer",n="integer",al12="numeric",
                         al13="numeric",al23="numeric",bal12="numeric",bal13="numeric",bal23="numeric",pbs="numeric")
likeCPP_code_HudsonFst<-"
double* q1 = new double[n[0]];
double* q2 = new double[n[0]];
double* q3 = new double[n[0]];

for(int i=0;i<n[0];i++){
   
   q1[i]=1-f1[i];
   q2[i]=1-f2[i];
   q3[i]=1-f3[i];

   al12[i] = (f1[i]-f2[i])*(f1[i]-f2[i]) - ((f1[i]*q1[i])/(N1[0]-1)) - ((f2[i]*q2[i])/(N2[0]-1));
   bal12[i] = f1[i]*q2[i] + f2[i]*q1[i];

   al13[i] = (f1[i]-f3[i])*(f1[i]-f3[i])* - ((f1[i]*q1[i])/(N1[0]-1)) - ((f3[i]*q3[i])/(N3[0]-1));
   bal13[i] = f1[i]*q3[i] + f3[i]*q1[i];

   al23[i] = (f2[i]-f3[i])*(f2[i]-f3[i]) - ((f2[i]*q2[i])/(N2[0]-1)) - ((f3[i]*q3[i])/(N3[0]-1));
   bal23[i] = f2[i]*q3[i] + f3[i]*q2[i];

  } 
for(int i=0;i<n[0];i++){
  pbs[i]= (-log(1-(al12[i]/bal12[i])) + -log(1-(al13[i]/bal13[i])) - -log(1-(al23[i]/bal23[i]))) / 2.0;
      
}
"

fns <- cfunction(list(pbsCalculator_HudsonFst=likeCPP_input_HudsonFst),
                 list(likeCPP_code_HudsonFst),
                 convention=".C", cxxargs="-O3", cppargs="-O3",language="C++")
pbsCalculator_forCpp_HudsonFst<-fns[["pbsCalculator_HudsonFst"]]











# 
# n<-length(right_pos)
# al12<-numeric(n)
# al13<-numeric(n)
# al23<-numeric(n)
# bal12<-numeric(n)
# bal13<-numeric(n)
# bal23<-numeric(n)
# pbs<-numeric(n)
# 
# system.time(cppPBS<-pbsCalculator_forCpp(fall[right_pos,"PEL"],fall[right_pos,"CEU"],fall[right_pos,"CHB"],nInd["PEL"],nInd["CEU"],nInd["CHB"],pos[right_pos],chr[right_pos],
#                              n=n,al12=al12,al13=al13,al23=al23,bal12=bal12,bal13=bal13,bal23=bal23,pbs=pbs))
# 
# plot(cppPBS[["pos"]],cppPBS[["pbs"]])
