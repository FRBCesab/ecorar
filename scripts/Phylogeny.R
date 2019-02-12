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
library("magrittr")

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

# Dropping names not in mammals ID
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

 
          
# Second version of the graph      
          rownames(data_DR)<-as.character(gsub(" ", "_", data_DR$Name))
          #order in the same order of the phylo
          data_DR<-data_DR[set_mammals$tip.label,] 

          data_DR$cols <- NA
          data_DR$cols[data_DR$DR_class=="D75R75"] <- "red"
          data_DR<-data_DR[rownames(data_DR) %in% mammalsPhy$tip.label,]
          
          rarety <-  data_DR$cols
          rarety <- as.numeric(as.factor(rarety))
          rarety[is.na(rarety)]<-0
          names(rarety)<-rownames(data_DR)
          
          # Change capita 
          names(nodesArc)<- dplyr::mutate_all(as.character(unique(data_DR$order)), funs=tolower)
          names(nodesArc) %<>% tolower
          
          firstup <- function(x) {
            substr(x, 1, 1) <- toupper(substr(x, 1, 1))
            x
          }
          names(nodesArc)<-firstup(names(nodesArc))
          
          
          # plotting PHYLOGENY TREE
          # plot(set_mammals,type = "fan",edge.color = "grey", show.tip.label = TRUE,tip.color="white", edge.width = 0.4)  # ,edge.lty= data_DR$ltystyle,edge.color = data_DR$cols)
          color.terminal.branches(set_mammals, rarety, breaks=4, cols=c("#A6A6A666","red"), edge.width=c(0.4), show.tip.label=TRUE, non.terminal.col= "#A6A6A666")
          tiplabels(pch = 18, col = data_DR$cols, cex = 0.4 ,offset=5)

          # plotting family labels/arcs
          offset <- rep(c(1.10,1.18,1.26),length(nodesArc)/2)
        
           # rep(c(1.08,1.16,1.24), length(nodesArc)/2)
          for(i in 1:length(nodesArc)){
            
            if(i %in% c(seq(1,length(nodesArc), by = 2))) laboffset <- 0.03
            if(i %in% c(seq(2,length(nodesArc), by = 2))) laboffset <- 0.03
            
            if(i %in% c(1,4,7,10,13,16,19,22)){  #If odd 
              arc.cladelabels(text= paste0(i,""), #paste0(names(nodesArc)[i])
                              node=nodesArc[i],
                              ln.offset=offset[i],
                              lab.offset=offset[i]+laboffset, 
                              cex = 0.6, 
                              colarc = "gray82",
                              lwd = 1, 
                              lty = 1, 
                              orientation = "curved",
                              mark.node = FALSE,col="gray82")}
              
              
              if(i %in% c(2,5,8,11,14,17,20,23)){  #If odd 
                arc.cladelabels(text= paste0(i,""), #paste0(names(nodesArc)[i])
                                node=nodesArc[i],
                                ln.offset=offset[i],
                                lab.offset=offset[i]+laboffset, 
                                cex = 0.6, 
                                colarc = "gray75",
                                lwd = 1, 
                                lty = 1, 
                                orientation = "curved",
                                mark.node = FALSE,col="gray75")}
                
                if(i %in% c(3,6,9,12,15,18,21)){  #If odd 
                  arc.cladelabels(text= paste0(i,""), #paste0(names(nodesArc)[i])
                                  node=nodesArc[i],
                                  ln.offset=offset[i],
                                  lab.offset=offset[i]+laboffset, 
                                  cex = 0.6, 
                                  colarc = "gray68",
                                  lwd = 1, 
                                  lty = 1, 
                                  orientation = "curved",
                                  mark.node = FALSE,col="gray68")}
           
             text(x=270,y=130-(i*10),labels=paste0(i,": ",names(nodesArc)[i]),cex=0.4)
          }
          

     

      
    # offset <- c(1.552727,1.645455,1.429091,1.181818,1.243636,1.738182,1.614545,1.305455,1.150909,1.398182,1.367273,1.707273,1.769091,1.583636,1.460000,
          #1.212727,1.274545,1.120000,1.505,1.645455, 1.336364,1.74,1.435)
   
          
#Compute Lambda to know if functional rare species are packaged           
         
              # LAMBDA IS NOT GOOD FOR BINARY INFO
            
           #Number of species per order 
           sp_order<- unlist(lapply(unique(data_DR$order), function(x){ 
                      l <- length(set_mammals$tip.label[which(data_DR$order == x)])
                      names(l) <- x
                      l} ))
            
            rarety <- data_DR$colpointD75R75
            rarety <- as.numeric(as.factor(rarety))
            rarety[is.na(rarety)]<-0
            names(rarety)<-rownames(data_DR)
            
            lambda<-phytools::phylosig(set_mammals, rarety, method="lambda", test=TRUE)   
          
#Compute FRITZ to know if functional rare species are packaged                 
            onction phylo.d() dans le package caper
            
            Faut utiliser au préalable comparative.data() pour formater les data 
            (traits, phylo etc....)
   


            #OTHER VERSION OF GRAPH
            
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
            
            
            if((i %% 2) == 0){  #If odd 
            arc.cladelabels(text="  ", #paste0(names(nodesArc)[i])
            node=nodesArc[i],
            ln.offset=offset[i],
            lab.offset=offset[i]+laboffset, 
            cex = 0.2, 
            colarc = "black",
            lwd = 3, 
            lty = 1, 
            orientation = "curved",
            mark.node = FALSE)
            
            }else{ 
            arc.cladelabels(text= "  ",#paste0(names(nodesArc)[i])
            node=nodesArc[i],
            ln.offset=offset[i],
            lab.offset=offset[i]+laboffset, 
            cex = 0.2, 
            colarc = "grey",
            lwd = 3, 
            lty = 1, 
            orientation = "curved",
            mark.node = FALSE)
            }
            }
            # TODO MAKE THIS AUTOMATICATLY Add figure of mammals order  
            require(png)
            rodentia <- readPNG(file.path(data_dir,"mammals","Images_order","rodentia.png"))              
            rasterImage(rodentia,5,180,80,240)
            
            lagomorpha <- readPNG(file.path(data_dir,"mammals","Images_order","Lagomorpha.png"))       
            rasterImage(lagomorpha,-250,70,-200,110,fg="grey")

     
