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
library(ape)
library("RColorBrewer")

#LOAD TRAITS MAPS AND DISTRIB----
reso="50km"
# Load traits and distrib 
load(file=file.path(data_dir,"mammals","mammalsID.RData"))
load(file=file.path(results_dir,"mammals","mammalstrait.RData"))
load(file=file.path(results_dir,"mammals",reso,"occ_mammals_list.RData"))
load(file=file.path(results_dir,"mammals",reso,"FR_mammals.RData"))

load(file=file.path(data_dir,"mammals","taxaInfo.RData"))


#----
mammalsID<-mammalsID[mammalsID$ID %in% rownames(mammalstrait),]
taxaInfo<-taxaInfo[taxaInfo$ID %in% mammalsID$ID,]


# Load Phylogeny
load(file=file.path(data_dir,"mammals","mammalsPhy.RData"))

# Dropping names not mammals ID
set_mammals <- ape::drop.tip(mammalsPhy,mammalsPhy$tip.label[!is.element(mammalsPhy$tip.label,as.character(gsub(" ", "_", mammalsID$Name)))])


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
      data_DR<-merge(data_DR,taxaInfo,by.x="row.names", by.y="ID")
      data_DR <- data_DR[,-c(11:13,15,16)]

      # First version of the graph
          rownames(data_DR)<-as.character(gsub(" ", "_", data_DR$Name))
          #order in the same order of the phylo
          data_DR<-data_DR[set_mammals$tip.label,] 
          colour <-  viridis::viridis(6,option="D")
          colour[6]<- "grey"
          data_DR$cols <- colour[as.factor(data_DR$DR_class)]
  
          data_DR<-data_DR[rownames(data_DR) %in% mammalsPhy$tip.label,]
          
                    #ltystyle <- c(2,2,1,2,1,2)
          #data_DR$ltystyle<- ltystyle[as.factor(data_DR$DR_class)]
     
          
          # Prepare family labels
          labelsArc <- as.character(unique(data_DR$order))
          labelsArc <- labelsArc[-which(labelsArc == "DERMOPTERA")]
          labelsArc <- labelsArc[-which(labelsArc == "TUBULIDENTATA")]
          labelsArc <- labelsArc[-which(labelsArc == "MICROBIOTHERIA") ]# the three order are absent from phylogeny but super weird like flying squirrel!!!
          
          nodesArc <- unlist(lapply(labelsArc, function(x){
            
            #x<-labelsArc[1]
            node <- phytools::findMRCA(set_mammals,as.character(set_mammals$tip.label[which(as.character(data_DR$order) == x)]), type = "node")
            names(node) <- x
            node
          }))
          nodesArc <- nodesArc[order(nodesArc, decreasing = FALSE)]
          
          # plotting PHYLOGENY TREE
          plot(set_mammals,type = "fan",edge.color = "grey", show.tip.label = TRUE,tip.color="white", edge.width = 0.4) # ,edge.lty= data_DR$ltystyle,edge.color = data_DR$cols)
          data_DR$colpointD75R75 <- data_DR$cols
          data_DR$colpointD75R75[ data_DR$colpointD75R75!="#7AD151FF"]<- NA
          
          data_DR$colpointD25R75 <- data_DR$cols
          data_DR$colpointD25R75[ data_DR$colpointD25R75!="#2A788EFF"]<- NA
          
          data_DR$colpointD25R25 <- data_DR$cols
          data_DR$colpointD25R25[ data_DR$colpointD25R25!="#414487FF"]<- NA
          
          tiplabels(pch = 19, col = data_DR$colpointD75R75, cex = 0.4 ,offset=5)
          tiplabels(pch = 19, col = data_DR$colpointD25R75, cex = 0.4 ,offset=10)
          tiplabels(pch = 19, col = data_DR$colpointD25R25, cex = 0.4 ,offset=15)
          
          # plotting family labels/arcs
          # offset <- c(seq(1.01,2.7, by = ((2.7-1.01)/23))[1:4],rep(seq(1.01,2.7, by = ((2.7-1.01)/23))[5:7], length(nodesArc)/3))
          # offset <- offset+0.1
          offset <- rep(c(1.20,1.15), length(nodesArc)/2)
          
        for(i in 1:length(nodesArc)){
            
            if(i %in% c(seq(1,length(nodesArc), by = 2))) laboffset <- 0.06
            if(i %in% c(seq(2,length(nodesArc), by = 2))) laboffset <- 0.06
            
            
            if((i %% 2) == 0){ 
                 arc.cladelabels(text="  ",
                            node=nodesArc[i],
                            ln.offset=offset[i],
                            lab.offset=offset[i]+laboffset, 
                            cex = 0.1, 
                            colarc = "black",
                            lwd = 3, 
                            lty = 1, 
                            orientation = "curved",
                            mark.node = FALSE)
                 
            }else{ 
                 arc.cladelabels(text="  ",
                                 node=nodesArc[i],
                                 ln.offset=offset[i],
                                 lab.offset=offset[i]+laboffset, 
                                 cex = 0.1, 
                                 colarc = "grey",
                                 lwd = 3, 
                                 lty = 1, 
                                 orientation = "curved",
                                 mark.node = FALSE)
            }
       }
          
          
          

                          
    
          
          
          
          
          
          data_DR$test<-NA
          data_DR$test[data_DR$DR_class=="D75R75"] <- 4
          data_DR$test[data_DR$DR_class=="D25R25"] <- 3
          data_DR$test[data_DR$DR_class=="D25R75"] <- 2
          data_DR$test[data_DR$DR_class=="D75R25"] <- 1
          data_DR$test[is.na(data_DR$test)] <- 0
          

          
          plotBranchbyTrait(set_mammals, x=data_DR$test, mode="edge", palette="rainbow", type = "fan")
          
          plotBranchbyTrait(tree, x, mode=c("edges","tips","nodes"), palette="rainbow", 
                            legend=TRUE, xlims=NULL, ...)
          
          
          
          
          
          
          
          
# Add figure  
          require(png)
          img<-readPNG(file=file.path(data_dir,"mammals",Images order, Afrosoricida.png))              
                        
                        
          #now open a plot window with coordinates
          plot(1:10,ty="n")
          #specify the position of the image through bottom-left and top-right coords
          rasterImage(img,2,2,4,4)          
          
          
          
          
          
          

          
          
        



#plot images of 1 representant of each family (clade)

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


                 
                 
draw.arc(50, 50, 100, deg2 = 1:10*10, col = "blue")





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











# Second version

#colour <-c("cadetblue1","orange","red2")
#data_DR[data_DR== "AVG"] <- NA
#data_DR[data_DR== "D75R25"] <- NA
#ata_DR[data_DR== "NA"] <- NA
#
data_DR<-merge(data_DR,mammalsID,by.x="row.names",by.y="ID")
rownames(data_DR)<-as.character(gsub(" ", "_", data_DR$Name))
#order in the same order of the phylo
data_DR<-data_DR[set_mammals$tip.label,] 

colour <-c("cadetblue1","orange","red2")
data_DR[data_DR== "AVG"] <- NA
data_DR[data_DR== "D75R25"] <- NA
data_DR[data_DR== "NA"] <- NA
data_DR$cols <- colour[as.factor(data_DR$DR_class)]
data_DR$cols[is.na(data_DR$cols)] <- "white"

plot(set_mammals,type = "fan",edge.color = "black", show.tip.label = FALSE, edge.width = 0.3)
tiplabels(pch = 15, col = data_DR$cols, cex = 1 )

