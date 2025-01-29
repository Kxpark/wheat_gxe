library(nasapower)
library(tidyverse)
library(EnvRtype)
library(reshape2)

locations<-data.frame('envs'=c('BD','BI','CH','ELS','PRO','EI','MCG'),
                      'lat'=c(35.1871,35.1871,34.1944,32.2672,33.1407,35.9303,31.3621),
                      'lon'=c(-102.0907,-102.0907,-99.5180,-96.7082,-96.0370,-101.9814,-97.4548))
years<-c(2019:2023)
env.i=c()
lat=c()
lon=c()
plant.date<-c()
harv.date<-c()
for (yr in years){
  
  env.i<-c(env.i,paste0(yr,locations$envs))
  lat<-c(lat,locations$lat)
  lon<-c(lon,locations$lon)
  yrhv<-yr-1
  plant.date<-as.Date(c(plant.date,rep(paste0(yrhv,"-10-30"),length(locations$envs))))
  harv.date<-as.Date(c(harv.date,rep(paste0(yr,"-06-30"),length(locations$envs))))
  
}

Data<-c()
for(envs in 1:length(lat)){

daily_test<-get_power(community = 'ag',
                      lonlat = c(lon[envs],lat[envs]),
                      pars=c("PRECTOTCORR","GWETTOP","TS","GWETROOT","T2M","T2M_MAX","T2M_MIN","WS2M","RH2M","T2MDEW","ALLSKY_SFC_LW_DWN","ALLSKY_SFC_SW_DWN"),
                      dates=c(plant.date[envs],harv.date[envs]),
                      temporal_api="daily")
colnames(daily_test)<-sub("PRECTOTCORR","PRECTOT",colnames(daily_test))
daily_test<-daily_test%>%mutate('daysFromStart'=1:nrow(daily_test))
daily_test<-daily_test%>%mutate("env"=env.i[envs])
daily_test<-as.data.frame(daily_test)
daily_test = processWTH(env.data = daily_test,Tbase1 = 0,Tbase2 = 45,Topt1 = 25,Topt2 = 28) ### set for wheat
daily_test<-daily_test%>%mutate("CumuPrecp"=cumsum(daily_test$PRECTOT),"CumuGDD"=cumsum(daily_test$GDD))
Data<-rbind(Data,daily_test)
}
write.csv(Data,file="ClimateData2019to2023.csv")
Soilparms<-read.delim("soilavg.txt")
Soilparms<-dcast(Soilparms,value.var = 'value',Env~Variable)
