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
###same order for rownames 
funk_mammals<-funk_mammals[match(map@data[,1], rownames(funk_mammals)),]

##PLOT maps
varmap <- names(funk_mammals)[-1]
lapply(varmap,function(i) map.Funk(data=funk_mammals,map=map,var=i,nlevels=5,plotpdf=TRUE,resultdir="mammals/50km",dalto=TRUE))
map.Funk(data=funk_mammals,map=map,var=varmap[3],nlevels=10,plotpdf=FALSE,resultdir="mammals",dalto=FALSE)




#----

#----

#BIRDS Rarity----

## LOAD birds values of functional rarity indices
load(file.path(results_dir,"birds","funk_birds.RData"))

###organised rownames 
funk_birds$cell<-as.numeric(as.character(funk_birds$cell))
rownames(funk_birds)<-funk_birds$cell

###Link MapGrid with value of indices
###same order for rownames 
funk_birds<-funk_birds[match(map@data[,1], rownames(funk_birds)),]

##PLOT maps
varmap <- names(funk_birds)[-1]
lapply(varmap,function(i) map.Funk(data=funk_birds,map=map,var=i,nlevels=10,plotpdf=TRUE,resultdir="birds/50km",dalto=FALSE))
map.Funk(data=funk_birds,map=map,var=varmap[3],nlevels=10,plotpdf=FALSE,resultdir="birds",dalto=FALSE)


#----

#FIND adjacent polygons in R (neighbors) ----
library(rgeos)
library(rgdal)
library(sp)
library(raster)
combined <- info %>%
  group_by(chuncks) %>%
  do(sf::st_intersection(., map))

a <- adjacent(map, ID_cell, pairs=F)

##Neighbour for map_mammals
Mat_neighbour<-gTouches(map, byid=T)
###Transform in 1/0
Mat_neighbour<-Mat_neighbour*1
save(Mat_neighbour,file=file.path(data_dir,"ReferenceGrid50km"))

neighbor<-mclapply(ID_cell[1:3],function(x) adjacent(map,x , directions=4, pairs=TRUE, target=NULL, sorted=FALSE, 
         include=FALSE, id=FALSE),mc.cores=3)


0.1)
