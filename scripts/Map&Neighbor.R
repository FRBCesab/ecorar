##############################################################################################################################
# MAPPING FUNCTION RARITY
#
# Author : Nicolas Mouquet, Nicolas Loiseau, 
# Date : 10-04-2018
# 
#R version 3.4.1 (2017-06-30) -- "Single Candle"
#Copyright (C) 2017 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#
#
#DATA from Thuiller W.
## 6 Apr 2018: spatial layers for global birds and .  are IUCN range maps, birds are breeding ranges only from bird life
## grid is 110km2 equal area with >70% water cells exluded on land, but oceanic islands with <70% water retained
###############################################################################################################################
rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NL")

library(parallel)
library(funrar)
library(ggplot2)

# Load spatial grid for plotting
map<-readOGR(file.path(data_dir,"ReferenceGrid50Km","gridLand50km.shp"))
#names of each cell
ID_cell<-map@data[,1]
#Resolution
reso="50km"

#mammals Rarity----

## LOAD birds values of functional rarity indices
load(file.path(results_dir,"mammals",reso,"funk_mammals.RData"))
load(file=file.path(results_dir,"mammals",reso,"FR_mammals.RData"))
load(file=file.path(results_dir,"mammals",reso,"sub_mammals.RData"))

###organised rownames 
funk_mammals$cell<-as.numeric(as.character(funk_mammals$cell))
rownames(funk_mammals)<-funk_mammals$cell

###Link MapGrid with value of indices
###delete polygons without any species information
##map_mammals<-map[map@data[,1] %in% funk_mammals$cell,]

##funk_mammals<-funk_mammals[funk_mammals$cell %in% map@data[,1],]

###same order for rownames 
funk_mammals<-funk_mammals[match(map@data[,1], rownames(funk_mammals)),]

##PLOT maps
varmap <- names(funk_mammals)[-1]
lapply(varmap,function(i) map.Funk(data=funk_mammals,map=map_mammals,var=i,nlevels=10,plotpdf=TRUE,resultdir="mammals",dalto=FALSE))
map.Funk(data=funk_mammals,map=map_mammals,var=varmap[3],nlevels=10,plotpdf=FALSE,resultdir="mammals",dalto=FALSE)
map.Funk(data=funk_mammals,map=map_mammals,var=varmap[6],nlevels=10,plotpdf=FALSE,resultdir="mammals",dalto=FALSE)

#----

#----

#BIRDS Rarity----

## LOAD birds values of functional rarity indices
load(file.path(results_dir,"birds","funk_birds.RData"))

###organised rownames 
funk_birds$poly<-as.numeric(as.character(funk_birds$poly))
rownames(funk_birds)<-funk_birds$poly

###Link MapGrid with value of indices
###delete polygons without any species information
map_birds<-priorities.shape[priorities.shape$HBWID %in% funk_birds$poly,]
funk_birds<-funk_birds[funk_birds$poly %in% map_birds$HBWID,]

###same order for rownames 
funk_birds<-funk_birds[match(map_birds$HBWID, rownames(funk_birds)),]

##PLOT maps

varmap <- names(funk_birds)[-1]
lapply(varmap,function(i) map.Funk(data=funk_birds,map=map_birds,var=i,nlevels=10,plotpdf=TRUE,resultdir="birds",dalto=FALSE))

#----

#FIND adjacent polygons in R (neighbors) ----
library(rgeos)
library(rgdal)
library(sp)

combined <- info %>%
  group_by(chuncks) %>%
  do(sf::st_intersection(., map))


##Neighbour for map_mammals
Mat_neighbour_mammals<-gTouches(map, byid=TRUE)
###Transform in 1/0
Mat_neighbour_mammals<-Mat_neighbour_mammals*1
save(Mat_neighbour_mammals,file=file.path(results_dir,"mammals","Mat_neighbour_mammals.RData"))
##Neighbour for map_birds
Mat_neighbour_birds<-gTouches(map_birds, byid=TRUE)
###Transform in 1/0
Mat_neighbour_birds<-Mat_neighbour_birds*1
save(Mat_neighbour_birds,file=file.path(results_dir,"birds","Mat_neighbour_birds.RData")) 
#----


#CHECK and VISUALISE neighbour ----     

plot(map_birds,lwd=0.1)

load(file.path(results_dir,"birds","Mat_neighbour_birds.RData"))

spot<-rownames(Mat_neighbour_birds)[sample(1:dim(Mat_neighbour_birds)[1],1)]
spot <- '19545'
spot_voisin<-Mat_neighbour_birds[spot,]
spot_voisin<-spot_voisin[spot_voisin>0]

plot(map_birds[map_birds$HBWID==spot,],col="red",add=T,lwd=0.1)
for (i in 1:length(spot_voisin)){
  plot(map_birds[map_birds$HBWID==names(spot_voisin)[i],],col="blue",add=T,lwd=0.1)
}

