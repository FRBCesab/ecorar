
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



#Grid ----

# Load spatial grid for plotting
grid<-readOGR(file.path(data_dir,"ReferenceGrid10Km","gridLand10km.B.shp"))
#names of each cell
ID_cell<-rownames(grid@data)

#Mammals----

#Occurence 
load(file=file.path(data_dir,"mamals/occ_mammals_sparseM.RData"))
#Check if each species have at least on occurence
occ_mammals <- occ_mammals[,colSums(occ_mammals)>0]

#ID 
load(file=file.path(data_dir,"mamals/MammalsID.RData"))
mammalsID<-data.frame(Mammals)
rm(Mammals)

mammalsID$checkname<-NA
for (i in 1:length(mammalsID$Name)){
  mammalsID$checkname[i] <- as.matrix(gnr_resolve(names = as.character(mammalsID$Name[i])))[1,3]
}
mammalsID$checkname[5128] <- as.character("Pseudoberylmys muongbangensis")

save(mammalsID,)




#


length(mammalsID$checkname%in%mamalstrait$checkname)

#Birds----

#Occurence 
load(file=file.path(data_dir,"mamals/occ_birds_sparseM.RData"))
#Check if each species have at least on occurence
occ_birds <- occ_birds[,colSums(occ_birds)>0]

#ID 
load(file=file.path(data_dir,"mamals/birdsID.RData"))
birdsID<-data.frame(birds)
rm(birds)

birdsID$checkname<-NA
for (i in 1:length(birdsID$Name)){
  birdsID$checkname[i] <- as.matrix(gnr_resolve(names = as.character(birdsID$Name[i])))[1,3]
}        





