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
## grid is 50km2 equal area with >70% water cells exluded on land, but oceanic islands with <70% water retained
###############################################################################################################################
rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NL")

library(parallel)
library(funrar)
library(ggplot2)
library(rgdal)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(latticeExtra)
library(devtools)
source_url("https://raw.githubusercontent.com/neocarto/R/master/colors/Palettes.R")
load("R/colors/colors.RData")

# Load spatial grid for plotting
map<-readOGR(file.path(data_dir,"ReferenceGrid50Km","gridLand50km.shp"))
#names of each cell
ID_cell<-map@data[,1]

#mammals Rarity----

## LOAD birds values of functional rarity indices
load(file.path(results_dir,"mammals","50km","funk_mammals.RData"))
load(file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))
load(file=file.path(results_dir,"mammals","50km","sub_mammals.RData"))
load(file=file.path(results_dir,"mammals","50km","SES_total_mammals.RData"))


###organised rownames 
funk_mammals$cell<-as.numeric(as.character(funk_mammals$cell))
rownames(funk_mammals)<-funk_mammals$cell

###Link MapGrid with value of indices
###same order for rownames 
funk_mammals<-funk_mammals[match(map@data[,1], rownames(funk_mammals)),]



##########################
#--- MAP MAMMALS

      # Join data with shapefile
      mapData <- merge(map, funk_mammals, by.x = 'ID', by.y = 'cell')
      mapDataNull <- merge(map, SES_total_mammals, by.x = 'ID', by.y = 'cell')
      
      # Import the world shapefile 
      World<-readOGR(file.path(data_dir,"ReferenceGrid50Km","Continents","GSHHS_i_L1.shp"))
      
      # Delete islands polygons (easy to read)
      World2 <- World[World@data$area > 20000,]
      
      # Change projection
      World2<-spTransform(World2,proj4string(mapData))
      
      # Plot Mammals
      quantiles_D75R75<-c(0,1,2,3,4,5,max(mapData@data["D75R75"]))
      mapData@data$cutD75R75 <- cut(mapData@data$D75R75, breaks = c(quantiles_D75R75),include.lowest = T,right=F) 
      
      quantiles_D25R25<-c(0,1,5,10,15,20,max(mapData@data["D25R25"]))
      mapData@data$cutD25R25 <- cut(mapData@data$D25R25, breaks = c(quantiles_D25R25),include.lowest = T,right=F) 
      
      quantiles_TD_sp<-c(0,1,10,20,50,100,max(mapData@data["TD_sp"]))
      mapData@data$cutTD_sp <- cut(mapData@data$TD_sp, breaks = c(quantiles_TD_sp),include.lowest = T,right=F) 
      
      quantiles_D75R75_Null <- c(round(min(na.omit(mapDataNull@data$D75R75)),1),-4,-1.96,1.96,4,round(max(na.omit(mapDataNull@data$D75R75)),1))
      mapDataNull@data$cutD75R75_Null  <- cut(mapDataNull@data$D75R75, breaks = c(quantiles_D75R75_Null),include.lowest = T,right=F) 
      
      pal_colNull <- c("#1874CDB3","#00CED1B3","#FFF68FB3","#FFA500B3","#FF0000B3")
      pal_col<- c("#EDEDEDB3", "#1874CDB3", "#00CED1B3", "#FFF68FB3", "#FFA500B3","#FF0000B3")
      
      map1 <-spplot(mapData["cutTD_sp"],col.regions = pal_D75R75,main = "TD_sp",
                    ## set the border color and width
                    col="transparent",
                    #col = pal[1], lwd = 0.01,
                    par.settings = list(axis.line=list(col="transparent")),
                    colorkey = list(height = 1, space = 'right',
                                    labels = list(at = seq(0.5, length(quantiles_TD_sp) -0.5),
                                                  labels = quantiles_TD_sp)),
                    contour = T) + layer(sp.polygons(World2, lwd = 0.6))
      
      
      map2 <-spplot(mapData["cutD75R75"],col.regions = pal_D75R75,main = "D25R25",
                    ## set the border color and width
                    col="transparent",
                    #col = pal[1], lwd = 0.01,
                    par.settings = list(axis.line=list(col="transparent")),
                    colorkey = list(height = 1, space = 'right',
                                    labels = list(at = seq(0.5, length(quantiles_D75R75) -0.5),
                                                  labels = quantiles_D75R75)),
                    contour = T) + layer(sp.polygons(World2, lwd = 0.6))
      
      map3 <-spplot(mapData["cutD25R25"],col.regions = pal_D25R25,main = "D25R25",
                    ## set the border color and width
                    col="transparent",
                    #col = pal[1], lwd = 0.01,
                    par.settings = list(axis.line=list(col="transparent")),
                    colorkey = list(height = 1, space = 'right',
                                    labels = list(at = seq(0.5, length(quantiles_D25R25) -0.5),
                                                  labels = quantiles_D25R25)),
                    contour = T) + layer(sp.polygons(World2, lwd = 0.6))
      
      map4 <-spplot(mapDataNull["cutD75R75_Null"],col.regions = pal_DataNull,main = "SES D75R75",
                    ## set the border color and width
                    col="transparent",
                    #col = pal[1], lwd = 0.01,
                    par.settings = list(axis.line=list(col="transparent")),
                    colorkey = list(height = 1, space = 'right',
                                    labels = list(at = seq(0.5, length(quantiles_D75R75_Null) -0.5),
                                                  labels = quantiles_D75R75_Null)),
                    contour = T) + layer(sp.polygons(World2, lwd = 0.6))
      
      pdf(file.path(results_dir,"mammals","50km",paste0("figs"),paste0("map","testAllmap",".pdf")))
      grid.arrange(map1,map2,map3,map4,nrow=4)
      dev.off()
      
      

#--- MAP BIRDS

## LOAD birds values of functional rarity indices
load(file.path(results_dir,"birds","50km","funk_birds.RData"))
load(file=file.path(results_dir,"birds","50km","SES_total_birds.RData"))
###organised rownames 
funk_birds$cell<-as.numeric(as.character(funk_birds$cell))
rownames(funk_birds)<-funk_birds$cell

###Link MapGrid with value of indices
###same order for rownames 
funk_birds<-funk_birds[match(map@data[,1], rownames(funk_birds)),]
##PLOT maps

###organised rownames 
funk_birds$cell<-as.numeric(as.character(funk_birds$cell))
rownames(funk_birds)<-funk_birds$cell

###Link MapGrid with value of indices
###same order for rownames 
funk_birds<-funk_birds[match(map@data[,1], rownames(funk_birds)),]

# Join data with shapefile
mapData <- merge(map, funk_birds, by.x = 'ID', by.y = 'cell')
mapDataNull <- merge(map, SES_total_birds, by.x = 'ID', by.y = 'cell')


# Plot Mammals
quantiles_D75R75<-c(0,1,2,5,10,20,max(mapData@data["D75R75"]))
mapData@data$cutD75R75 <- cut(mapData@data$D75R75, breaks = c(quantiles_D75R75),include.lowest = T,right=F) 

quantiles_D25R25<-c(0,1,5,20,50,70,max(mapData@data["D25R25"]))
mapData@data$cutD25R25 <- cut(mapData@data$D25R25, breaks = c(quantiles_D25R25),include.lowest = T,right=F) 

quantiles_TD_sp<-c(0,1,10,100,300,500,max(mapData@data["TD_sp"]))
mapData@data$cutTD_sp <- cut(mapData@data$TD_sp, breaks = c(quantiles_TD_sp),include.lowest = T,right=F) 

quantiles_D75R75_Null <- c(round(min(na.omit(mapDataNull@data$D75R75)),1),-4,-1.96,1.96,round(max(na.omit(mapDataNull@data$D75R75)),1))
mapDataNull@data$cutD75R75_Null  <- cut(mapDataNull@data$D75R75, breaks = c(quantiles_D75R75_Null),include.lowest = T,right=F) 

pal_DataNull <- c("#1874CDB3","#00CED1B3","#FFF68FB3","#FF0000B3")
pal_D75R75<- c("#EDEDEDB3", "#1874CDB3", "#00CED1B3", "#FFF68FB3", "#FFA500B3","#FF0000B3")

map1 <-spplot(mapData["cutTD_sp"],col.regions = pal_D75R75,
              ## set the border color and width
              col="transparent",
              #col = pal[1], lwd = 0.01,
              par.settings = list(axis.line=list(col="transparent")),
              colorkey = list(height = 1, space = 'right',
                              labels = list(at = seq(0.5, length(quantiles_TD_sp) -0.5),
                                            labels = quantiles_TD_sp)),
              contour = T) + layer(sp.polygons(World2, lwd = 0.6))


map2 <-spplot(mapData["cutD75R75"],col.regions = pal_D75R75,
              ## set the border color and width
              col="transparent",
              #col = pal[1], lwd = 0.01,
              par.settings = list(axis.line=list(col="transparent")),
              colorkey = list(height = 1, space = 'right',
                              labels = list(at = seq(0.5, length(quantiles_D75R75) -0.5),
                                            labels = quantiles_D75R75)),
              contour = T) + layer(sp.polygons(World2, lwd = 0.6))

map3 <-spplot(mapData["cutD25R25"],col.regions = pal_D25R25,
              ## set the border color and width
              col="transparent",
              #col = pal[1], lwd = 0.01,
              par.settings = list(axis.line=list(col="transparent")),
              colorkey = list(height = 1, space = 'right',
                              labels = list(at = seq(0.5, length(quantiles_D25R25) -0.5),
                                            labels = quantiles_D25R25)),
              contour = T) + layer(sp.polygons(World2, lwd = 0.6))

map4 <-spplot(mapDataNull["cutD75R75_Null"],col.regions = pal_DataNull,
              ## set the border color and width
              col="transparent",
              #col = pal[1], lwd = 0.01,
              par.settings = list(axis.line=list(col="transparent")),
              colorkey = list(height = 1, space = 'right',
                              labels = list(at = seq(0.5, length(quantiles_D75R75_Null) -0.5),
                                            labels = quantiles_D75R75_Null)),
              contour = T) + layer(sp.polygons(World2, lwd = 0.6))

pdf(file.path(results_dir,"birds","50km",paste0("figs"),paste0("map","testAllmap",".pdf")))
grid.arrange(map1,map2,map3,map4,nrow=4)
dev.off()










##PLOT maps OLD VERSION 
varmap <- names(funk_mammals)[-1]
lapply(varmap,function(i) map.Funk(data=funk_mammals,map=map,var=i,nlevels=5,plotpdf=TRUE,resultdir="mammals/50km",dalto=TRUE))
map.Funk(data=funk_mammals,map=map,var=varmap[3],nlevels=5,plotpdf=FALSE,resultdir="mammals",dalto=FALSE)

varmap <- names(funk_birds)[-1]
lapply(varmap,function(i) map.Funk(data=funk_birds,map=map,var=i,nlevels=10,plotpdf=TRUE,resultdir="birds/50km",dalto=FALSE))
map.Funk(data=funk_birds,map=map,var=varmap[3],nlevels=10,plotpdf=FALSE,resultdir="birds",dalto=FALSE)


