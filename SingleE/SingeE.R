library(BGLR)
library(tidyverse)
library(corrplot)
library(doParallel)
pheno<-read.csv("/scratch/user/kxparker/GS/BM/phenoReady.csv") #pheno

pheno<-pheno[order(pheno$Environment),]
GM <- read.delim("/scratch/user/kxparker/GS/BM/Final Test/GRM.txt") #GRM
reps<-20
Pop<-"TXBreeding19to23"
envs<-unique(pheno$Environment)

PredictionSingle<-function(x,pheno,GM,envs,reps){
    P2<-pheno%>%filter(Environment==envs[x])
    GM2<-as.matrix(GM[P2$Genotype,P2$Genotype])
    ETM<-list(GM=list(K=GM2,model="RKHS")) #you can change the model here
    CORFOLD<-matrix(nrow=reps,ncol=1,dimnames=list(c(1:reps),envs[x]))
 for(i in 1:reps){
     corFld<-matrix(nrow=5,ncol=1,dimnames = list(c(1:5), envs[x]))
     P3<-P2
     Fold<-sample(1:5,nrow(P3),replace=T)
       for(z in 1:5){
    P3<-P2$BLUE_yield
    P3[which(Fold==z)]<-NA
    fm<-BGLR(y=P3,ETA = ETM,nIter = 5000,burnIn = 2000)
  
    corFld[z,1]<-cor(P2[which(Fold==z),'BLUE_yield'],fm$yHat[which(Fold==z)] ,use = "complete.obs")
   }#end of fold
    CORFOLD[i,1]<-mean(corFld[,1],na.rm=T)
     
 } #end of rep  
    return(CORFOLD)
    
}#end of function




cl<-makeCluster(24)
registerDoParallel(cl)
DataFinal<-foreach(x=1:length(envs),.combine='cbind',.packages=c('tidyverse','BGLR'))%dopar%{PredictionSingle(x,pheno,GM,envs,reps)}
write.table(DataFinal,file="Corfile.SE_solo_RKHS.txt")
stopCluster(cl)
