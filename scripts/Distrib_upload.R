
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
# map<-readOGR(file.path(data_dir,"ReferenceGrid10Km","gridLand10km.B.shp"))
#names of each cell
# ID_cell<-rownames(map@data)

#Mammals----

load(file.path(results_dir,"mammals/mammalstrait.RData"))

    #Occurence 
    load(file=file.path(data_dir,"mammals/occ_mammals_sparseM.RData"))

    #ID 
    load(file=file.path(results_dir,"mammals/MammalsID.RData"))

    #Commun ID for mammalsID/occ_mammals/traitmammals ---
    mammalsID<-mammalsID[mammalsID$checkname_clean %in% mammalstrait$checkname_clean,]
    mammalstrait<-merge(mammalsID,mammalstrait,by="checkname_clean")
    rownames(mammalstrait)<-mammalstrait$ID
    mammalstrait<-mammalstrait[,-c(2:7,23,24)]
    
    occ_mammals <- occ_mammals[,colnames(occ_mammals)  %in% mammalsID$ID]
    save(occ_mammals, file=file.path(results_dir,"mammals/occ_mammals.RData"))
    #Transform Matrix occ in list to reduce objectsize
    occ_mammals_list <- lapply(occ_mammals, function(.col){sp <- row.names(occ_mammals)[.col == 1]})

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





