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

###organised rownames 
funk_mammals$cell<-as.numeric(as.character(funk_mammals$cell))
rownames(funk_mammals)<-funk_mammals$cell

###Link MapGrid with value of indices
###same order for rownames 
funk_mammals<-funk_mammals[match(map@data[,1], rownames(funk_mammals)),]

##PLOT maps
varmap <- names(funk_mammals)[-1]
lapply(varmap,function(i) map.Funk(data=funk_mammals,map=map,var=i,nlevels=5,plotpdf=TRUE,resultdir="mammals/50km",dalto=TRUE))
map.Funk(data=funk_mammals,map=map,var=varmap[3],nlevels=5,plotpdf=FALSE,resultdir="mammals",dalto=FALSE)


# SORTIR LES POLYGONES TROP PETIT
# MODIF COULEUR


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




##########################
# Import the world shapefile 
World<-readOGR(file.path(data_dir,"ReferenceGrid50Km","Continents","GSHHS_i_L1.shp"))

# Delete islands polygons (easy to read)
World2 <- World[World@data$area > 20000,]

# Change projection
World2<-spTransform(World2,proj4string(mapData))


#--- MAP MAMMAMLS

      # Join data with shapefile
      mapData <- merge(map, funk_mammals, by.x = 'ID', by.y = 'cell')


    # Compute breaks with quantiles and assign color
          
       # Number of classes for the age data
    
          #D25R25
          
          no_classes_D25R25  <- 11
          quantiles_D25R25 <- quantile(funk_mammals["D25R25"], 
                                       probs = seq(0, 1,
                                                   length.out = no_classes_D25R25),na.rm=T)
          pal_D25R25 <- rev(brewer.pal(no_classes_D25R25, "Spectral"))
          pal_D25R25<- pal_D25R25[c(1:5,7:11)]
          pal_D25R25<-c("white",pal_D25R25) #gray92
          #displaypal(pal_D25R25)

          #D75R75
          
          no_classes_D75R75 <- 6
          quantiles_D75R75 <- quantile(funk_mammals["D75R75"], 
                                       probs = seq(0, 1,
                                                   length.out = no_classes_D75R75),na.rm=T)
          quantiles_D75R75 <- c(0,1,2,3,4,13)
          pal_D75R75<-pal_D25R25[c(1,7,8,9,10,11)]

          #TD_sp
          
          no_classes_TD_sp <- 11
          quantiles_TD_sp <- quantile(funk_mammals["TD_sp"], 
                                      probs = seq(0, 1,
                                                  length.out = no_classes_TD_sp),na.rm=T)
          quantiles_TD_sp[2] <- 1
          pal_TD_sp <- pal_D25R25

map1 <- spplot(mapData["D75R75"],
              col.regions = pal_D75R75,main = "D75R75",
              ## define the points where the colors change
              at = quantiles_D75R75,
              ## set the border color and width
              col="transparent",
              #col = pal[1], lwd = 0.01,
              par.settings = list(axis.line=list(col="transparent")),
              ## adjust the legend
              #colorkey =
              #  list(space = 'right',
              #       height = 0.3,
              #       labels = list(
              #         at = quantiles_D75R75,
              #         labels = signif(quantiles_D75R75, 1),
              #         rot = 0))
                )+ layer(sp.polygons(World2, lwd = 0.6))

map2 <- spplot(mapData["TD_sp"],
               col.regions = pal_TD_sp,main = "TD_sp",
               ## define the points where the colors change
               at = quantiles_TD_sp,
               ## set the border color and width
               col="transparent",
               par.settings = list(axis.line=list(col="transparent")),
               ## adjust the legend
               #colorkey =
               # list(space = 'right',
               #      height = 0.3,
               #      labels = list(
               #        at = quantiles_TD_sp,
               #        labels = signif(quantiles_TD_sp, 1),
               #        rot = 0))
                 )+ layer(sp.polygons(World2, lwd = 0.6))

map3 <- spplot(mapData["D25R25"],
               col.regions = pal_D25R25,main = "D25R25",
               ## define the points where the colors change
               at = quantiles_D25R25,
               ## set the border color and width
               col="transparent",
               par.settings = list(axis.line=list(col="transparent")),
               ## adjust the legend
               #colorkey =
               # list(space = 'right',
               #      height = 0.3,
               #      labels = list(
               #        at = quantiles_D25R25,
               #        labels = signif(quantiles_D25R25, 1),
               #        rot = 0))
                 )+ layer(sp.polygons(World2, lwd = 0.6))

pdf(file.path(results_dir,"mammals","50km",paste0("figs"),paste0("map","testAllmap",".pdf")))
grid.arrange(map2,map1,map3,nrow=3)
dev.off()






#--- MAP BIRDS

# Join data with shapefile
mapData <- merge(map, funk_birds, by.x = 'ID', by.y = 'cell')

# Compute breaks with quantiles and assign color

# Number of classes for the age data

#D25R25

no_classes_D25R25  <- 11
quantiles_D25R25 <- quantile(funk_birds["D25R25"], 
                             probs = seq(0, 1,
                                         length.out = no_classes_D25R25),na.rm=T)
quantiles_D25R25[2] <- 1
pal_D25R25 <- rev(brewer.pal(no_classes_D25R25, "Spectral"))
pal_D25R25<- pal_D25R25[c(1:5,7:11)]
pal_D25R25<-c("white",pal_D25R25) #gray92
#displaypal(pal_D25R25)

#D75R75

no_classes_D75R75 <- 6
quantiles_D75R75 <- quantile(funk_birds["D75R75"], 
                             probs = seq(0, 1,
                                         length.out = no_classes_D75R75),na.rm=T)
quantiles_D75R75 <- c(0,1,4,8,12,18)
pal_D75R75<-pal_D25R25[c(1,7,8,9,10,11)]

#TD_sp

no_classes_TD_sp <- 11
quantiles_TD_sp <- quantile(funk_birds["TD_sp"], 
                            probs = seq(0, 1,
                                        length.out = no_classes_TD_sp),na.rm=T)
quantiles_TD_sp[1]<-0
quantiles_TD_sp[2] <- 1
pal_TD_sp <- pal_D25R25

map1 <- spplot(mapData["D75R75"],
               col.regions = pal_D75R75,main = "D75R75",
               ## define the points where the colors change
               at = quantiles_D75R75,
               ## set the border color and width
               col="transparent",
               #col = pal[1], lwd = 0.01,
               par.settings = list(axis.line=list(col="transparent")),
               ## adjust the legend
               #colorkey =
               #   list(space = 'right',
               #    height = 0.3,
               #      labels = list(
               #         at = quantiles_D75R75,
               #        labels = signif(quantiles_D75R75, 3),
               #         rot = 0))
                 )+ layer(sp.polygons(World2, lwd = 0.6))

map2 <- spplot(mapData["TD_sp"],
               col.regions = pal_TD_sp,main = "TD_sp",
               ## define the points where the colors change
               at = quantiles_TD_sp,
               ## set the border color and width
               col="transparent",
               par.settings = list(axis.line=list(col="transparent")),
               ## adjust the legend
               #colorkey =
               #  list(space = 'right',
               #       height = 0.3,
               #       labels = list(
               #         at = quantiles_TD_sp,
               #         labels = signif(quantiles_TD_sp, 1),
               #         rot = 0))
                 )+ layer(sp.polygons(World2, lwd = 0.6))

map3 <- spplot(mapData["D25R25"],
               col.regions = pal_D25R25,main = "D25R25",
               ## define the points where the colors change
               at = quantiles_D25R25,
               ## set the border color and width
               col="transparent",
               par.settings = list(axis.line=list(col="transparent")),
               ## adjust the legend
               #colorkey =
               #  list(space = 'right',
               #       height = 0.3,
               #       labels = list(
               #         at = quantiles_D25R25,
               #         labels = signif(quantiles_D25R25, 1),
               #         rot = 0))
                 )+ layer(sp.polygons(World2, lwd = 0.6))

pdf(file.path(results_dir,"birds","50km",paste0("figs"),paste0("map","testAllmap",".pdf")))
grid.arrange(map2,map1,map3,nrow=3)
dev.off()

