library(BGLR)
library(tidyverse)
library(doParallel)
library(matrixcalc)

pheno<-read.csv("/scratch/user/kxparker/GS/BM/phenoReady.csv") #pheno
pheno<-pheno[order(pheno$Environment),]
gm <- read.delim("/scratch/user/kxparker/GS/BM/Final Test/GRM.txt") #GRM
gm<-as.matrix(gm)

#Select a W matrix
#wenv<-as.matrix(read.delim("/scratch/user/kxparker/GS/BM/wrmReady.txt",sep=' '))
#wenv<-as.matrix(read.delim("/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W0matrix.txt",sep='\t')) #W0
#wenv<-as.matrix(read.delim("/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W1matrix.txt",sep='\t')) #W1
#wenv<-as.matrix(read.delim("/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W2matrix.txt",sep='\t')) #W2
wenv<-as.matrix(read.delim("/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W3matrix.txt",sep='\t')) #W3
#wenv<-as.matrix(read.delim("/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W4matrix.txt",sep='\t')) #W4
#wenv<-as.matrix(read.delim("/scratch/user/kxparker/GS/BM/Final Test/Wmatrixs/W5matrix.txt",sep='\t')) #W5


colnames(wenv)<-rownames(wenv)
wenv<-wenv[unique(pheno$Environment),unique(pheno$Environment)]





P1<-data.frame()
for (i in 1:nrow(wenv)){
  P1<-rbind(P1,data.frame('Environment'=rownames(wenv)[i],'Genotype'=colnames(gm)))
}

Zy<-diag(1,nrow(wenv),nrow(wenv))

Kronc<-kronecker(Zy,gm) # all by all matrix
#Kronc<-kronecker(wenv,gm)

pheno2<-pheno %>% mutate('Combo'=paste0(Environment,Genotype))%>%select(Combo,BLUE_yield)
P1<-P1%>%mutate('Combo'=paste0(Environment,Genotype))
P1<-merge(P1,pheno2,by='Combo',all.x=T)
P1<-P1%>%filter(!duplicated(Combo))
P2<-P1[-which(is.na(P1$BLUE_yield)),]
ENV<-as.factor(P2$Environment)
Ze<-model.matrix(~ENV-1)
Kronc<-Kronc[-which(is.na(P1$BLUE_yield)),-which(is.na(P1$BLUE_yield))]

ETM<-list(Gkronc=list(X=Kronc,model='BRR'),
           Envs=list(X=Ze,model='FIXED'))
#ETM<-list(Gkronc1=list(X=Kronc,model='BRR'))
#ETM<-list(Gkronc1=list(K=Kronc,model='RKHS'),
#            Envs=list(X=Ze,model='FIXED'))


P3<-P2%>%select(BLUE_yield)
ENV<-as.factor(P2$Environment)
envs<-unique(P2$Environment)


reps<-10
Pop<-"TXBreeding19to23"



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
LOOgsmatrixSingleRep<-function(x,P2,ETM,envs,ENV){
  Cormatrix<-matrix(nrow=1,ncol = 1,dimnames = list(1,envs[x])) 
  P3<-P2$BLUE_yield
  P3[which(ENV==envs[x])]<-NA
  fm<-BGLR(y=P3,ETA = ETM,groups = ENV,nIter = 5000,burnIn = 2000)  
  Cormatrix[1,1]<-cor(P2[which(ENV==envs[x]),'BLUE_yield'],fm$yHat[which(ENV==envs[x])], use = "complete.obs")
  write.table(Cormatrix,file=paste0("Cormatrix.data.",x,".single.Kronctest.txt"))
  return(Cormatrix)
} 



cl<-makeCluster(12)
registerDoParallel(cl)

#DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{LOOgsmatrixRKHS(x,P2,ETM,envs,ENV)}
#write.table(DataFinal,file="Cormatrix.GxE_W3_RKHS_fixed.LOO.txt")
DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{LOOgsmatrix(x,P2,ETM,envs,ENV)}
write.table(DataFinal,file="Cormatrix.GxE_fixed_BRR_fixed.LOO.txt")

#DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{PredictionRKHS(x,P2,ETM,envs,ENV)}
#write.table(DataFinal,file="Corfile.GxE_fixed_RKHS_fixed.txt")
#DataFinal<-foreach(x=1:reps,.combine='rbind',.packages=c('tidyverse','BGLR'))%dopar%{Prediction(x,P2,ETM,envs,ENV)}
#write.table(DataFinal,file="Corfile.GxE_fixed_BRR_fixed.txt")
stopCluster(cl)


