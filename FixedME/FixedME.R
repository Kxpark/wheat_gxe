library(BGLR)
library(tidyverse)
library(corrplot)
library(doParallel)
library(matrixcalc)
pheno<-read.csv("/scratch/user/kxparker/GS/BM/phenoReady.csv") #pheno

pheno<-pheno[order(pheno$Environment),]
GM <- as.matrix(read.delim("/scratch/user/kxparker/GS/BM/Final Test/GRM.txt")) #GRM



Zg<-model.matrix(~pheno$Genotype-1)
envs<-unique(pheno$Environment)
ENV<-as.factor(pheno$Environment)
Ze<-model.matrix(~ENV-1)
colnames(Ze)<-sub("ENV","",colnames(Ze))

GM<-Zg%*%GM%*%t(Zg)

reps<-10
cl<-makeCluster(20)
registerDoParallel(cl)

#for(Wmatrix in c(FALSE,"wrm","W0","W1","W2","W3","W4","W5")){
for(Wmatrix in c("W2","W3","W4","W5")){
P2<-pheno

if(Wmatrix!=FALSE){
    #Select a W matrix
    Wm<-c("/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/wrmReady2.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W0matrix.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W1matrix.txt",
    "/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W2matrix.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W3matrix.txt","/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W4matrix.txt",
    "/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W5matrix.txt")
    names(Wm)<-c("wrm","W0","W1","W2","W3","W4","W5")
    wenv<-as.matrix(read.delim(Wm[Wmatrix],sep='\t'))
    filename<-paste0("Corfile.ME_",Wmatrix,"_BRR_LOO.txt")
    colnames(wenv)<-rownames(wenv)
    wenv<-wenv[envs,envs]
    
    Wmx<-Ze%*%wenv%*%t(Ze)
    
    ETM<-list(GM1=list(X=GM,model="BRR"),
            WX=list(X=Wmx,model="BRR"))
    
    
}else{
    filename<-paste0("Corfile.ME_fixed_BRR_LOO.txt")
    ETM<-list(GM1=list(X=GM,model="BRR"),
            WX=list(X=Ze,model="FIXED"))
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

#DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{PredictionRKHS(x,P2,ETM,envs,ENV)}
#DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{Prediction(x,P2,ETM,envs,ENV)}
#DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{LOOgsmatrixRKHS(x,P2,ETM,envs,ENV)}
DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{LOOgsmatrix(x,P2,ETM,envs,ENV)}

write.table(DataFinal,file=paste0(filename))
}
stopCluster(cl)
