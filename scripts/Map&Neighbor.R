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
library(rgdal)
# Load spatial grid for plotting
map<-readOGR(file.path(data_dir,"ReferenceGrid50Km","gridLand50km.shp"))
#names of each cell
ID_cell<-map@data[,1]

#mammals Rarity----

## LOAD birds values of functional rarity indices
load(file.path(results_dir,"mammals","50km","funk_mammals.RData"))
load(file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))
load(file=file.path(results_dir,"mammals","50km","sub_mammals.RData"))

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
load(file.path(results_dir,"birds","50km","funk_birds.RData"))

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







#--- New function to map functional rarity
#column <- colnames(funk_mammals)[c(2:7)]
#log_funk_mammals <- do.call(cbind,lapply(column, function(i) log10(funk_mammals[i])+1))
#log_funk_mammals <- cbind(funk_mammals[,1],log_funk_mammals)
#colnames(log_funk_mammals)[1] <- "cell"
#log_funk_mammals[log_funk_mammals==-Inf]<-0

## Join data with shapefile
mapData <- merge(map, funk_mammals, by.x = 'ID', by.y = 'cell')
## Import the raster file 

## Number of classes for the age data
no_classes <- 6

## Colors using the magma palette

pal <-  magma(n = 200,
              ## Limit the first and last colors
              begin = 0.1,
              end = 0.8)

#pal<-pal[c(1,5,10,15,20,25,30,35)]
pal<-pal[c(1,25,50,75,100,125,150)]
pal<-rev(pal)
pal[1]<-"gray92"
pal[2]<-"gold" 
pal[4]<-pal[3]
pal[3]<-"orange2"
pal[7]<-"black"

displaypal(pal)
## Compute breaks with quantiles
quantiles <- quantile(mapData$D75R75, 
                      probs = seq(0, 1,
                                  length.out = no_classes + 1),na.rm=T)

if(var=="D75R75") quantiles[1:8] <- c(0,1,2,3,4,6,8,10)

## Not needed here, but a better alternative to original code (a "for"
## loop)
## labels <- paste(signif(quantiles[-7], 4),
##                 signif(quantiles[-1], 4),
##                 sep = '-')

## Create a trellis object using sp::spplot to display the
## SpatialPolygons object.spplot
pdf(file.path(results_dir,resultdir,"50km",paste0("figs"),paste0("map","testD75R75",".pdf")))
spplot(mapData["D75R75"],
              col.regions = pal,
              ## define the points where the colors change
              at = quantiles,
              ## set the border color and width
              col = "#F2F2F202", lwd = 0.02,
              par.settings = list(axis.line=list(col="transparent")),
              ## adjust the legend
              colorkey =
                list(space = 'bottom',
                     height = 0.7,
                     labels = list(
                       at = quantiles,
                       labels = signif(quantiles, 3),
                       rot = 30)
                )
)

dev.off()






## Join data with shapefile
mapData <- merge(map, log_funk_mammals, by.x = 'ID', by.y = 'cell')
## Import the raster file 

## Number of classes for the age data
no_classes <- 40

## Colors using the magma palette

pal <-  magma(n = no_classes + 1 ,
              ## Limit the first and last colors
              begin = 0.1,
              end = 0.8)

pal<-pal[c(3,9,12,15,18,21,24,27,30,33,36)]
pal<-rev(pal)
pal[1]<-"gray92"
pal[2]<-"goldenrod1" 
pal[4]<-pal[3]
pal[3]<-"orange2"
pal[11]<-"black"
displaypal(pal)
## Compute breaks with quantiles
quantiles <- quantile(mapData$D75R75, 
                      probs = seq(0, 1,
                                  length.out = 10 + 1),na.rm=T)



quantiles[1:11]<-c(0,0.21,0.42,0.63,0.84,1.05,1.26,1.47,1.68,1.89,2.1)





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
