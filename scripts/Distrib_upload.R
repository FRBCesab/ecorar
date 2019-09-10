
rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NL")
library(funrar)
library(moments)
library(ggplot2)
library(parallel)
library(ade4)
library(dplyr)
library(gridExtra)
library(cluster)
library(rgdal)
library(Matrix)
library(taxize)

if (!require("devtools")) install.packages("devtools")
#Choose resolution----
reso="50km"

#Mammals----

load(file.path(results_dir,"mammals",reso,"mammalstrait.RData"))
#ID 
load(file=file.path(data_dir,"mammals","mammalsID.RData"))

    #Commun ID for mammalsID/occ_mammals/traitmammals ---
    mammalsID<-mammalsID[mammalsID$checkname_clean %in% mammalstrait$checkname_clean,]
    mammalstrait<-merge(mammalsID,mammalstrait,by="checkname_clean")
    rownames(mammalstrait)<-mammalstrait$ID
    mammalstrait<-mammalstrait[,-c(1:7,23,24)]
    save(mammalstrait, file=file.path(results_dir,"mammals",reso,"mammalstrait.RData"))
    
    
    load(file=file.path(data_dir,"mammals",reso,"occ_mammals_list.RData"))
    occ_mammals_list<-mammalsPresence
    rm(mammalsPresence)
    #
    #Work with a list, delete species that are not in mammalsID
    for (i in 1: length(occ_mammals_list)){
      occ_mammals_list[[i]]<- occ_mammals_list[[i]][occ_mammals_list[[i]]%in%mammalsID$ID]
      print(paste0('i',i))
    }
    save(occ_mammals_list,file=file.path(results_dir,"mammals",reso,"occ_mammals_list.RData"))
    
    
    
    
#Birds----

#Occurence 
load(file=file.path(data_dir,"birds",reso,"birdsPresence160119.RData"))
load(file=file.path(data_dir,"birds",reso,"birdsID.RData"))

birdstrait<-read.csv2(file.path(data_dir,"birds",reso,"BirdFuncDat.csv"))
birdstrait<-birdstrait[birdstrait$Scientific %in% birdsID[,2],]

birdstrait<-merge(birdsID,birdstrait,by.x="Name",by.y="Scientific") 
rownames(birdstrait)<-birdstrait$ID
birdstrait<-birdstrait[,c(11:20,25:32,36,37)]
birdstrait$BodyMass.Value<-as.numeric(as.character(birdstrait$BodyMass.Value))
save(birdstrait,file=file.path(results_dir,"birds",reso,"birdstrait.RData"))

birdsID <- birdsID[birdsID[,1] %in% rownames(birdstrait),]
save(birdsID,file=file.path(results_dir,"birds",reso,"birdsID.RData"))

#ID 
occ_birds_list<-birdsPresence
rm(birdsPresence)
#
#Work with a list, delete species that are not in birdsID
for (i in 1: length(occ_birds_list)){
  occ_birds_list[[i]]<- occ_birds_list[[i]][occ_birds_list[[i]]%in%birdsID[,1]]
  print(paste0('i',i))
}

names(occ_birds_list)<-names(occ_mammals_list)
save(occ_birds_list,file=file.path(results_dir,"birds",reso,"occ_birds_list.RData"))






