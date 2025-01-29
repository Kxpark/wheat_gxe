library(BGLR)
library(tidyverse)
pheno<-read.csv("/scratch/user/kxparker/GS/BM/phenoReady.csv") #pheno
pheno<-pheno[order(pheno$Environment),]
GM <- read.delim("/scratch/user/kxparker/GS/BM/Final Test/GRM.txt") #GRM
GM<-GM[pheno$Genotype,pheno$Genotype]
reps<-20
Pop<-"TXBreeding19to23"
envs<-unique(pheno$Environment)

## MultiEnv with GxE##
#set up#

ENV<-as.factor(pheno$Environment)
Ze<-model.matrix(~ENV-1)
colnames(Ze)<-sub("ENV","",colnames(Ze))
P2<-pheno
#ETM<-list(Envs=list(X=Ze,model="FIXED"),
#          GM=list(X=GM,model="BRR"))
ETM<-list(GM2=list(X=GM,model="BRR"))

#start pairwise Env testing, unknown env tests#
Cormatrix<-matrix(nrow=reps,ncol = length(envs),dimnames = list(c(1:reps),envs))
 for(z in 1:reps){ 
  for (i in 1:length(envs)){
    P3<-P2$BLUE_yield
    P3[which(ENV==envs[i])]<-NA
    fm<-BGLR(y=P3,ETA = ETM,groups = ENV,nIter = 5000,burnIn = 2000)  
    Cormatrix[z,i]<-cor(P2[which(ENV==envs[i]),5],fm$yHat[which(ENV==envs[i])], use = "complete.obs")
  }#end of matrix
  write.table(Cormatrix,file="Cormatrix.ME_solo.LOO.txt")
 }
 
unlink("*.dat") 
 
#start pairwise Env testing, unknown env tests#
#Cormatrix<-matrix(nrow=length(envs),ncol = length(envs),dimnames = list(envs,envs))
#  for (i in 1:length(envs)){
#    P3<-P2$BLUE_yield
#    P3[which(ENV!=envs[i])]<-NA
#    fm<-BGLR(y=P3,ETA = ETM,groups = ENV,nIter = 5000,burnIn = 2000)
#    for (M in 1:length(envs)){
#      Cormatrix[i,M]<-cor(P2[which(ENV==envs[M]),5],fm$yHat[which(ENV==envs[M])], use = "complete.obs")
#      write.table(Cormatrix,file="Cormatrix.fixed.txt")
#    }#end of rep
#  }#end of matrix
