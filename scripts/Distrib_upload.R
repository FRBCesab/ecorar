
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

#map ----

# Load spatial grid for plotting
map<-readOGR(file.path(data_dir,"ReferenceGrid10Km","gridLand10km.B.shp"))
#names of each cell
ID_cell<-rownames(map@data)

#Mammals----

load(file.path(results_dir,"mammals/mammalstrait.RData"))


    #Occurence 
    load(file=file.path(data_dir,"mammals/occ_mammals_sparseM.RData"))


    #ID 
    load(file=file.path(data_dir,"mammals/MammalsID.RData"))
    mammalsID<-data.frame(Mammals)
    rm(Mammals)# renames the files 
    
    #Check names ---
    #mammalsID$checkname<-NA
    #for (i in 1:length(mammalsID$Name)){
    #mammalsID$checkname[i] <- as.matrix(gnr_resolve(names = as.character(mammalsID$Name[i])))[1,3]
    #print(paste0("i",i))
    #}
    #mammalsID$checkname[5128] <- as.character("Pseudoberylmys muongbangensis")
    
    #mammalsID$Name<-mammalsID$checkname
    #mammalsID<-mammalsID[,-3]
    #colnames(mammalsID)[2]<-"checkname"
    #save(mammalsID, file=file.path(results_dir,"mammals/mammalsID.RData"))
    
    load(file=file.path(results_dir,"mammals/mammalsID.RData"))
    #Commun ID for mammalsID/occ_mammals/traitmammals ---
    mammalsID<-mammalsID[mammalsID$checkname %in% mammalstrait$checkname,]
    mammalstrait<-merge(mammalsID,mammalstrait,by="checkname")
    rownames(mammalstrait)<-mammalstrait$ID
    mammalstrait<-mammalstrait[,-c(2,3)]
    
    occ_mammals <- occ_mammals[,colnames(occ_mammals)  %in% mammalsID$ID]
    #Check if each species have at least on occurence
    occ_mammals <- occ_mammals[,colSums(occ_mammals)>0]

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





