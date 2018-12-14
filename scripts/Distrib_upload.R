
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
load(file=file.path(data_dir,"mammals",reso,"MammalsID.RData"))

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
load(file=file.path(data_dir,"mammals/occ_birds_sparseM.RData"))
#Check if each species have at least on occurence
occ_birds <- occ_birds[,colSums(occ_birds)>0]

#ID 
load(file=file.path(data_dir,"mammals/birdsID.RData"))
birdsID<-data.frame(birds)
rm(birds)

birdsID$checkname<-NA
for (i in 1:length(birdsID$Name)){
  birdsID$checkname[i] <- as.matrix(gnr_resolve(names = as.character(birdsID$Name[i])))[1,3]
}        





