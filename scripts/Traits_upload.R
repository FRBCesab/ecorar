############
#DATA UPLOAD
#Metadata for the traits : http://www.esapubs.org/archive/ecol/E095/178/metadata.php
############

rm(list=ls(all=TRUE)) 

library(parallel)
library(ggplot2)
library(taxize)

source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NL")

#MAMALS Traits ----

mammalstrait_all<- read.csv(file.path(data_dir,"mammals","MamFuncDat.csv"),sep=";") 

mammalstrait_all$Diet.Inv <- as.numeric(mammalstrait_all$Diet.Inv)
mammalstrait_all$Diet.Vend <- as.numeric(mammalstrait_all$Diet.Vend)
mammalstrait_all$Diet.Vect <- as.numeric(mammalstrait_all$Diet.Vect)
mammalstrait_all$Diet.Vfish <- as.numeric(mammalstrait_all$Diet.Vfish)
mammalstrait_all$Diet.Vunk <- as.numeric(mammalstrait_all$Diet.Vunk)
mammalstrait_all$Diet.Scav <- as.numeric(mammalstrait_all$Diet.Scav)
mammalstrait_all$Diet.Fruit <- as.numeric(mammalstrait_all$Diet.Fruit)
mammalstrait_all$Diet.Nect <- as.numeric(mammalstrait_all$Diet.Nect)
mammalstrait_all$Diet.Seed <- as.numeric(mammalstrait_all$Diet.Seed)
mammalstrait_all$Diet.PlantO <- as.numeric(mammalstrait_all$Diet.PlantO)
mammalstrait_all$ForStrat.Value<- as.factor(mammalstrait_all$ForStrat.Value)
mammalstrait_all$Activity.Nocturnal<- as.factor(mammalstrait_all$Activity.Nocturnal)
mammalstrait_all$Activity.Crepuscular<- as.factor(mammalstrait_all$Activity.Crepuscular)
mammalstrait_all$Activity.Diurnal<- as.factor(mammalstrait_all$Activity.Diurnal)
mammalstrait_all$BodyMass.Value<- as.numeric(mammalstrait_all$BodyMass.Value)

interest.var.mammals <- c("MSW3_ID","Scientific","Diet.Inv","Diet.Vend","Diet.Vect","Diet.Vfish","Diet.Vunk","Diet.Scav","Diet.Fruit","Diet.Nect","Diet.Seed",
                         "Diet.PlantO","ForStrat.Value","Activity.Nocturnal","Activity.Crepuscular","Activity.Diurnal","BodyMass.Value")

mammalstrait <- mammalstrait_all[,interest.var.mammals]

#Checknames = Verified names of spaces
mammalstrait$checkname<-NA
for (i in 1:length(mammalstrait$Scientific)){
  mammalstrait$checkname[i] <- as.matrix(gnr_resolve(names = as.character(mammalstrait$Scientific[i])))[1,3]
  print(paste0("i",i))
}
mammalstrait$Scientific<-mammalstrait$checkname
mammalstrait<-mammalstrait[,-18]
colnames(mammalstrait)[2]<-"checkname"
save(mammalstrait, file=file.path(results_dir,"mammals/mammalstrait.RData"))

