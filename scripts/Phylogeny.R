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
          data_DR<-data_DR[rownames(data_DR) %in% mammalsPhy$tip.label,]
          #order in the same order of the phylo
          data_DR<-data_DR[set_mammals$tip.label,] 
          
          data_DR$cols <- NA
          data_DR$cols[data_DR$DR_class=="D75R75"] <- "red"
    
          
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
       
          # Change capita 
          names(nodesArc)<- dplyr::mutate_all(as.character(unique(data_DR$order)), funs=tolower)
          names(nodesArc) %<>% tolower
          
          firstup <- function(x) {
            substr(x, 1, 1) <- toupper(substr(x, 1, 1))
            x
          }
          names(nodesArc)<-firstup(names(nodesArc))
         
           # Add column binary for Functional Rarity: Yes/no
          rarety <-  data_DR$cols 
          rarety <- as.numeric(as.factor(rarety))
          rarety[is.na(rarety)]<-0
          names(rarety)<-rownames(data_DR)
          
          set_mammals <- ape::drop.tip(mammalsPhy,mammalsPhy$tip.label[!is.element(mammalsPhy$tip.label,as.character(gsub(" ", "_", rownames(data_DR))))])
          
          # plotting PHYLOGENY TREE
          color.terminal.branches(set_mammals, rarety, breaks=4, cols=c("#A6A6A666","red"), edge.width=0.4, show.tip.label=TRUE, non.terminal.col= "#A6A6A666")
          tiplabels(pch = 18, col = data_DR$cols, cex = 0.4 ,offset=5)

          # plotting family labels/arcs
          offset <- rep(c(1.10,1.18,1.26),length(nodesArc)/2)
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
           
           #Add the names of Order
           text(x=270,y=130-(i*10),labels=paste0(i,": ",names(nodesArc)[i]),cex=0.4)
          }
          
#---
#Compute Lambda to know if functional rare species are packaged           

      #Compute FRITZ to know if functional rare species are packaged        
            data_DR$rarety <- data_DR$cols
            data_DR$rarety <- as.numeric(as.factor(data_DR$rarety))
            data_DR$rarety[is.na(data_DR$rarety)]<-0

    
          #Should be done with the 100 trees of mammals et 100 trees of birds
            data_DR<-data_DR
            data_DR$species<-rownames(data_DR)
            
D.phylogeny <- function(ids,proc,data_DR,taxa,permut) {
              #proc <- 2
              #taxa <- mammals
              #data_DR <- data_DR
              #permut <- 10
              #ids <- 1:2
  
mclapply(ids,function(id) { 
  if (taxa == "mammals") tree<-read.tree(file=file.path(data_dir,"mammals","alltrees", paste0("FritzTree_mammals_updateCarnivora2012DEF",id,".tre"))) 
  
  if (taxa == "birds")   tree<-read.tree(file=file.path(data_dir,"birds","alltrees", paste0("BirdzillaHackett1_",i,".tre"))) 
          set_tree<- ape::drop.tip(tree,tree$tip.label[!is.element(tree$tip.label,as.character(gsub(" ", "_", rownames(data_DR))))])
          set_tree$node.label <- NULL
        
            #collapse or resolve multichotomies in phylogenetic trees TODO check that is mean exactely because ned it
            set_tree<-di2multi(set_tree)
            #Compute D and statistic
            FR_PhyloD <- comparative.data(set_tree, data_DR,"species",na.omit=FALSE)
            FR_PhyloD <- phylo.d(FR_PhyloD, binvar=rarety,permut=1000)
            
            #The estimated D value
            estimated_D <- FR_PhyloD$DEstimate
            #A p value,giving the result of testing whether D is significantly different from one
            Pval1 <- FR_PhyloD$Pval1
            #A p value, giving the result of testing whether D is significantly different from zero
            Pval0 <- FR_PhyloD$Pval0
            
            },mc.cores= proc)
}
