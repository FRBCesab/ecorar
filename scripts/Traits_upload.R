############
#DATA UPLOAD
#Metadata for the traits : http://www.esapubs.org/archive/ecol/E095/178/metadata.php
############

rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NL")

library(parallel)
library(ggplot2)
library(taxize)


############
load(file=file.path(data_dir,"mammals/MammalsID.RData"))
mammalsID<-data.frame(Mammals)
rm(Mammals)


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


mammalsID$clean <- paste(word(mammalsID$Name,1),word(mammalsID$Name,2),sep=" ")
mammalstrait$clean <- paste(word(mammalstrait$Scientific,1),word(mammalstrait$Scientific,2),sep=" ")

#clean name
mammalsID$clean <- stringCleaning(mammalsID$clean)
mammalstrait$clean <- stringCleaning(mammalstrait$clean)

#Find duplicated species

sp_dup <- mammalsID[duplicated2(mammalsID[, "clean"]), ]#8 species are dupblicated due to subspecies (no functional traits for them, delete)

#Choose to delete spp

sp_dup<-sp_dup[grep("ssp.",sp_dup$Name, ignore.case = TRUE),]
mammalsID <- mammalsID[!rownames(mammalsID)%in%rownames(sp_dup),]


#Check on data base name
#mammalstrait$checkname<-NA
#for (i in 1:length(mammalstrait$clean)){
#  mammalstrait$checkname[i] <- as.matrix(gnr_resolve(names = as.character(mammalstrait$clean[i])))[1,3]
#  print(paste0("i",i))
#}

#mammalsID$checkname<-NA
#for (i in 5120:length(mammalsID$Name)){
#mammalsID$checkname[i] <- as.matrix(gnr_resolve(names = as.character(mammalsID$Name[i])))[1,3]
#print(paste0("i",i))
#}
#mammalsID$checkname[5119] <- as.character("pseudoberylmys muongbangensis")

#clean name
mammalsID$checkname_clean <- stringCleaning(mammalsID$checkname)
mammalstrait$checkname_clean <- stringCleaning(mammalstrait$checkname)


save(mammalstrait, file=file.path(results_dir,"mammals/mammalstrait.RData"))
save(mammalsID, file=file.path(results_dir,"mammals/mammalsID.RData"))