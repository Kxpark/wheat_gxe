#Working Sheet

setwd("C:/Users/cool_/OneDrive - Texas A&M AgriLife/Thesis/BM/R")
Data<-read.csv("ClimateFrom2010to2023.csv")[,-1]
#or
Data<-read.csv("ClimateData2019to2023.csv")[,-1]
soilparams<-read.delim("soilavg.txt")

library(tidyverse)
library(corrplot)
library(reshape2)
soilparams<-dcast(soilparams,value.var = 'value',Env~Variable)

D1<-as.matrix(Data[,-c(1:7,21)])

cor<-cor(D1)
corrplot(cor,order='hclust')
pc<-prcomp(D1)
pc$rotation<-pc$rotation*-1
biplot(pc)

biplot(pc)


D2<-Data%>%filter(daysFromStart>63)
stages<-c("Jan","Feb","March","April","May","June")
interval<-c(64,94,122,154,184,214)

EC10<-W_matrix(env.data = D2,env.id = 'env',var.id = c("FRUE","T2M_RANGE","PETP","PRECTOT","ETP","VPD","WS2M","SPV"),scale = T,center = T, by.interval = T,time.window = interval,
               names.window = stages,statistic = 'mean')
EC11<-cbind(EC10,"irr"=rep(0,nrow(EC10)))

EC11[grep("BI",rownames(EC11)),'irr']<-1
EC11[grep("EI",rownames(EC11)),'irr']<-.75
W0<-env_kernel(env.data = EC11,gaussian = T)[[2]]
write.table(W0,file="W0matrix.txt",sep='\t')
W1<-env_kernel(env.data = EC11,gaussian = F)[[2]]
write.table(W1,file="W1matrix.txt",sep='\t')

stages<-c("Nov","Dec","Jan","Feb","March","April","May","June")
interval<-c(0,33,64,94,122,154,184,214)

EC10<-W_matrix(env.data = Data,env.id = 'env',var.id = c("FRUE","T2M_RANGE","PETP","PRECTOT","ETP","VPD","WS2M","SPV","TS","CumuPrecp","CumuGDD","GWETROOT"),scale=F,center=F, by.interval = T,time.window = interval,
               names.window = stages,statistic = 'mean')
EC11<-cbind(EC10,"irr"=rep(0,nrow(EC10)))

EC11[grep("BI",rownames(EC11)),'irr']<-1
EC11[grep("EI",rownames(EC11)),'irr']<-.75

E10<-mutate("clay"=NA,"irr"=NA,"OC"=NA,"sand"=NA,"silt"=NA)


W2<-env_kernel(env.data = EC11,gaussian = F)[[2]]
write.table(W2,file="W2matrix.txt",sep='\t')
W3<-env_kernel(env.data = EC11,gaussian = T)[[2]]
write.table(W3,file="W3matrix.txt",sep='\t')

EC11<-cbind(EC10,"irr"=rep(0,nrow(EC10)),"clay"=rep(NA,nrow(EC10)),"OC"=rep(NA,nrow(EC10)),"sand"=rep(NA,nrow(EC10)),"silt"=rep(NA,nrow(EC10)))

for(env in rownames(soilparams)){
  EC11[grep(env,rownames(EC11)),'clay']<-soilparams[env,'clay']
  EC11[grep(env,rownames(EC11)),'irr']<-soilparams[env,'irr']
  EC11[grep(env,rownames(EC11)),'OC']<-soilparams[env,'OC']
  EC11[grep(env,rownames(EC11)),'sand']<-soilparams[env,'sand']
  EC11[grep(env,rownames(EC11)),'silt']<-soilparams[env,'silt']
  
  }

W4<-env_kernel(env.data = EC11,gaussian = F)[[2]]
write.table(W4,file="W4matrix.txt",sep='\t')
W5<-env_kernel(env.data = EC11,gaussian = T)[[2]]
write.table(W5,file="W5matrix.txt",sep='\t')




D2<-Data[which(Data$daysFromStart<154),]
stages<-c("Nov","Dec","Jan","Feb","March")
interval<-c(0,33,64,94,122)

EC10<-W_matrix(env.data = Data,env.id = 'env',var.id = c("FRUE","T2M_RANGE","PETP","PRECTOT","ETP","VPD","WS2M","SPV","TS","CumuPrecp","CumuGDD","GWETROOT"),scale=F,center=F, by.interval = T,time.window = interval,
               names.window = stages,statistic = 'mean')
EC11<-cbind(EC10,"irr"=rep(0,nrow(EC10)),"clay"=rep(NA,nrow(EC10)),"OC"=rep(NA,nrow(EC10)),"sand"=rep(NA,nrow(EC10)),"silt"=rep(NA,nrow(EC10)))


for(env in rownames(soilparams)){
  EC11[grep(env,rownames(EC11)),'clay']<-soilparams[env,'clay']
  EC11[grep(env,rownames(EC11)),'irr']<-soilparams[env,'irr']
  EC11[grep(env,rownames(EC11)),'OC']<-soilparams[env,'OC']
  EC11[grep(env,rownames(EC11)),'sand']<-soilparams[env,'sand']
  EC11[grep(env,rownames(EC11)),'silt']<-soilparams[env,'silt']
  
}

rownames(EC11)<-sub("2020","20",rownames(EC11))
rownames(EC11)<-sub("202","2",rownames(EC11))
rownames(EC11)<-sub("201","1",rownames(EC11))

W6<-env_kernel(env.data = EC11,gaussian = T)[[2]]
write.table(W6,file="W6matrix.txt",sep='\t')
W7<-env_kernel(env.data = EC11,gaussian = F)[[2]]
write.table(W7,file="W7matrix.txt",sep='\t')
