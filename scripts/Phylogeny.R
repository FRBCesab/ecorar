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
library("RColorBrewer")

#LOAD TRAITS MAPS AND DISTRIB----
reso="50km"
# Load traits and distrib 
load(file=file.path(data_dir,"mammals","mammalsID.RData"))
load(file=file.path(results_dir,"mammals","mammalstrait.RData"))
load(file=file.path(results_dir,"mammals",reso,"occ_mammals_list.RData"))
load(file=file.path(results_dir,"mammals",reso,"FR_mammals.RData"))
#----
mammalsID<-mammalsID[mammalsID$ID %in% rownames(mammalstrait),]

# Load Phylogeny
load(file=file.path(data_dir,"mammals","mammalsPhy.RData"))

#Add Information on Rarity
mammalsPhy[[mammalsPhy$tip.label %in% gsub(" ", "_", mammalsID$Name)]]


# Dropping names not mammals ID
set_mammals <- drop.tip(mammalsPhy,mammalsPhy$tip.label[!is.element(mammalsPhy$tip.label,as.character(gsub(" ", "_", mammalsID$Name)))])

# Create Class DR 

data_DR<-FR_mammals$FR
data_DR$DR_class="NA"

QD75 <- FR_mammals$Q$Q75_D
QD25 <- FR_mammals$Q$Q25_D
QR75 <- FR_mammals$Q$Q75_R
QR25 <- FR_mammals$Q$Q25_R

data_DR$DR_class[(data_DR$Din<QD25) & (data_DR$Rin<QR25)]="D25R25"
data_DR$DR_class[(data_DR$Din>QD75) & (data_DR$Rin>QR75)]="D75R75"
data_DR$DR_class[(data_DR$Din<QD25) & (data_DR$Rin>QR75)]="D25R75"
data_DR$DR_class[(data_DR$Din>QD75) & (data_DR$Rin<QR25)]="D75R25"
data_DR$DR_class[(((data_DR$Din>QD25) & (data_DR$Din<QD75)) & ((data_DR$Rin>QR25) & (data_DR$Rin<QR75)))]="AVG"

data_DR$InvRin=1-data_DR$Rin

# Create Class DR 
data_DR$order="NA"




data_DR<-merge(data_DR,mammalsID,by.x="row.names",by.y="ID")
rownames(data_DR)<-as.character(gsub(" ", "_", data_DR$Name))
data_DR<-data_DR[set_mammals$tip.label,] 

data_DR[data_DR== "AVG"] <- NA
data_DR[data_DR== "D75R25"] <- NA
data_DR[data_DR== "NA"] <- NA


colour <- rev(brewer.pal(3, "Spectral"))


data_DR$cols <- colour[as.factor(data_DR$DR_class)]
data_DR$cols[is.na(data_DR$cols)]<-"white"
plot(set_mammals,type = "fan",tip.color = data_DR$cols, cex = 0.1)

dend2 <- color_labels(dend, k = 3)

#plot images of 1 representant of each family (clade)
order <- 

  usedTree <- reorder(set_mammals, order = "cladewise")
# library(png)
w<-50
for(i in 1:length(spp)){
  # i = 13
  if(i %in% c(seq(1,length(spp), by = 2))) arrl <- 180
  if(i %in% c(seq(2,length(spp), by = 2))) arrl <- 180
  
  xy<-add.arrow(obj,spp[i],col="transparent",arrl=arrl,lwd=3,hedl=0.1)
  img<-load.image(file.path(pathphoto_png,paste(grep(sppp[i],list.files(pathphoto_png), value = TRUE)[1],sep="")))
  img <- as.raster(img)
  img[img=='#FFFFFF']=NA  # pour forcer la transparence
  asp<-dim(img)[1]/dim(img)[2]
  rasterImage(img,xy$x0-w/2,xy$y0-w/2*asp,xy$x0+w/2,xy$y0+w/2*asp)
  add.arrow(obj,spp[i],col="lightblue",arrl=arrl-22,lwd=2,hedl=0.05)
  rm(img)
}

# Prepare family labels
labelsArc <- as.character(unique(myFishesSerf$fam))
labelsUnq <- labelsArc[which(labelsArc == "Muraenidae")]
labelsArc <- labelsArc[-which(labelsArc == "Muraenidae")]
nodesArc <- unlist(lapply(labelsArc, function(x){
  
  node <- findMRCA(set100Serf[[100]],as.character(myFishesSerf$tip_name[which(as.character(myFishesSerf$fam) == x)]), type = "node")
  
  names(node) <- x
  node
}))

nodesArc <- nodesArc[order(nodesArc, decreasing = FALSE)]

# plotting family labels/arcs
offset <- c(seq(1.01,2.7, by = ((2.7-1.01)/17))[1:4],rep(seq(1.01,2.7, by = ((2.7-1.01)/17))[5:7], length(nodesArc)/3))
# offset <- c(sample(seq(1.1,1.8, by = ((2.3-1.1)/17))))
for(i in 1:length(nodesArc)){
  
  if(i %in% c(seq(1,length(nodesArc), by = 2))) laboffset <- 0.06
  if(i %in% c(seq(2,length(nodesArc), by = 2))) laboffset <- 0.06
  
  
  arc.cladelabels(text=paste0(names(nodesArc)[i]),
                  # "\n",
                  # round(mean(myFishesSerf$Mean_tot[which(as.character(myFishesSerf$fam) == names(nodesArc)[i])]),2),
                  # "(",
                  # "±",
                  # round(sd(myFishesSerf$Mean_tot[which(as.character(myFishesSerf$fam) == names(nodesArc)[i])]),2),
                  # ")"),
                  node=nodesArc[i],
                  ln.offset=offset[i],
                  lab.offset=offset[i]+laboffset, #+laboffset
                  cex = 0.7, 
                  col = sample(colors()[grep("gold",colors())], 1),
                  lwd = 1.3, 
                  orientation = "curved",
                  mark.node = FALSE)
}



arc.cladelabels(text=labelsUnq,node=which(obj$tree$tip.label=="Muraenidae_Gymnothorax_javanicus"),
                orientation="curved",ln.offset=1.25,lab.offset=1.32, cex = 0.7,col = sample(colors()[grep("gold",colors())], 1), mark.node = FALSE)

