library(BGLR)
library(tidyverse)
library(doParallel)
library(matrixcalc)

pheno<-read.csv("/scratch/user/kxparker/GS/BM/phenoReady.csv") #pheno
pheno<-pheno[order(pheno$Environment),]
gm <- read.delim("/scratch/user/kxparker/GS/BM/Final Test/GRM.txt") #GRM
gm<-as.matrix(gm)
envs<-unique(pheno$Environment)
LOO<-TRUE

if(LOO){
reps<-10
ending<-"LOO"
}else{
    reps=20
    ending<-"Cor"
}

#Wmatrix<-"W5" #FALSE "wrm","W0","W1","W2","W3","W4","W5"
cl<-makeCluster(12)
registerDoParallel(cl)
#c("wrm","W0","W1","W2","W3","W4","W5")
for(Wmatrix in c(FALSE,"W4","W5")){

if(Wmatrix!=FALSE){
    #Select a W matrix
    Wm<-c("/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/wrmReady2.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W0matrix.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W1matrix.txt",
    "/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W2matrix.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W3matrix.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W4matrix.txt",
    "/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W5matrix.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W6matrix.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W7matrix.txt")
    names(Wm)<-c("wrm","W0","W1","W2","W3","W4","W5","W6","W7")
    wenv<-as.matrix(read.delim(Wm[Wmatrix],sep='\t'))
    filename<-paste0("Corfile.GxE_",Wmatrix,"_handmard_nofix_",ending,".txt")
    colnames(wenv)<-rownames(wenv)
    wenv<-wenv[envs,envs]
}else{
    wenv<-diag(1,length(unique(pheno$Environment)),length(unique(pheno$Environment)))
    filename<-paste0("Corfile.GxE_fixed_handmard_nofix_",ending,".txt")
}

ENV<-as.factor(pheno$Environment)
Ze<-model.matrix(~ENV-1)
Zg<-model.matrix(~pheno$Genotype-1)


Geno<-Zg%*%gm%*%t(Zg)
Wenv<-Ze%*%wenv%*%t(Ze)

WxG<-hadamard.prod(Wenv,Geno)


#ETM<-list(Gkronc1=list(K=WxG,model='RKHS'),
#            Envs=list(X=Ze,model='FIXED'))
ETM<-list(Gkronc1=list(K=WxG,model='RKHS'))

P2<-pheno%>%select(BLUE_yield,Environment,Genotype)



PairwisePrediction<-function(x,P2,ETM,envs,ENV){
  Cormatrix<-matrix(nrow=1,ncol = length(envs),dimnames = list(envs[x],envs))
  P3<-P2$BLUE_yield
  P3[which(ENV!=envs[x])]<-NA
  fm<-BGLR(y=P3,ETA = ETM,groups = ENV,nIter = 5000,burnIn = 2000)
  for (M in 1:length(envs)){
    Cormatrix[1,M]<-cor(P2[which(ENV==envs[M]),'BLUE_yield'],fm$yHat[which(ENV==envs[M])], use = "complete.obs")}#end of rep
  
  return(Cormatrix)
}#end function




Prediction<-function(x,P2,ETM,envs,ENV){
  corFld<-matrix(nrow=5,ncol=length(envs),dimnames = list(c(1:5), envs))
  Fold<-sample(1:5,nrow(P2),replace = T)
  CorFEnv<-matrix(nrow=1,ncol=length(envs),dimnames = list(x, envs))
  for(z in 1:5){
    P3<-P2$BLUE_yield
    P3[which(Fold==z)]<-NA
    fm<-BGLR(y=P3,ETA = ETM,groups = ENV,nIter = 5000,burnIn = 2000)
    for(U in 1:length(envs)){
      corFld[z,U]<-cor(P2[which(which(Fold==z) %in% which(P2$Environment==envs[U])),'BLUE_yield'],fm$yHat[which(which(Fold==z) %in% which(P2$Environment==envs[U]))],use = "complete.obs")
    }#cor
    
  }#end of fold
  for(U in 1:length(envs)){
    CorFEnv[1,U]<-mean(corFld[,U],na.rm = T)}#end of means
  
  return(CorFEnv)
}

PredictionRKHS<-function(x,P2,ETM,envs,ENV){
  corFld<-matrix(nrow=5,ncol=length(envs),dimnames = list(c(1:5), envs))
  Fold<-sample(1:5,nrow(P2),replace = T)
  CorFEnv<-matrix(nrow=1,ncol=length(envs),dimnames = list(x, envs))
  for(z in 1:5){
    P3<-P2$BLUE_yield
    P3[which(Fold==z)]<-NA
    fm<-BGLR(y=P3,ETA = ETM,nIter = 5000,burnIn = 2000)
    for(U in 1:length(envs)){
      corFld[z,U]<-cor(P2[which(which(Fold==z) %in% which(P2$Environment==envs[U])),'BLUE_yield'],fm$yHat[which(which(Fold==z) %in% which(P2$Environment==envs[U]))],use = "complete.obs")
    }#cor
    
  }#end of fold
  for(U in 1:length(envs)){
    CorFEnv[1,U]<-mean(corFld[,U],na.rm = T)}#end of means
  
  return(CorFEnv)
}
LOOgsmatrixRKHS<-function(x,P2,ETM,envs,ENV){
  Cormatrix<-matrix(nrow=1,ncol = length(envs),dimnames = list(x,envs)) 
  for (i in 1:length(envs)){
    P3<-P2$BLUE_yield
    P3[which(ENV==envs[i])]<-NA
    fm<-BGLR(y=P3,ETA = ETM,nIter = 5000,burnIn = 2000)  
    Cormatrix[1,i]<-cor(P2[which(ENV==envs[i]),'BLUE_yield'],fm$yHat[which(ENV==envs[i])], use = "complete.obs")
  }#end of rep
  return(Cormatrix)
} 
LOOgsmatrix<-function(x,P2,ETM,envs,ENV){
  Cormatrix<-matrix(nrow=1,ncol = length(envs),dimnames = list(x,envs)) 
  for (i in 1:length(envs)){
    P3<-P2$BLUE_yield
    P3[which(ENV==envs[i])]<-NA
    fm<-BGLR(y=P3,ETA = ETM,groups = ENV,nIter = 5000,burnIn = 2000)  
    Cormatrix[1,i]<-cor(P2[which(ENV==envs[i]),'BLUE_yield'],fm$yHat[which(ENV==envs[i])], use = "complete.obs")
  }#end of rep
  return(Cormatrix)
} 

LOOgsmatrixSingleRep<-function(x,P2,ETM,envs,ENV){
  Cormatrix<-matrix(nrow=1,ncol = 1,dimnames = list(1,envs[x])) 
  P3<-P2$BLUE_yield
  P3[which(ENV==envs[x])]<-NA
  fm<-BGLR(y=P3,ETA = ETM,groups = ENV,nIter = 5000,burnIn = 2000)  
  Cormatrix[1,1]<-cor(P2[which(ENV==envs[x]),'BLUE_yield'],fm$yHat[which(ENV==envs[x])], use = "complete.obs")
  write.table(Cormatrix,file=paste0("Cormatrix.data.",x,".single.Kronctest.txt"))
  return(Cormatrix)
} 



if(LOO){
DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{LOOgsmatrixRKHS(x,P2,ETM,envs,ENV)}
}else{
DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{PredictionRKHS(x,P2,ETM,envs,ENV)}
}
write.table(DataFinal,file=paste0(filename))
}

stopCluster(cl)


