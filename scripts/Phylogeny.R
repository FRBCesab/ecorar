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
library(RColorBrewer)
library(magrittr)
library(ggplot2)
library(grid)
library(picante)

#LOAD TRAITS MAPS AND DISTRIB----

    # Load traits and distrib 

        #Mammals
        load(file=file.path(data_dir,"mammals","mammalsID.RData"))
        load(file=file.path(results_dir,"mammals","mammalstrait.RData"))
        load(file=file.path(results_dir,"mammals","50km","occ_mammals_list.RData"))
        load(file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))
        load(file=file.path(data_dir,"mammals","taxaInfo_mammals.RData"))
        mammalsID<-mammalsID[mammalsID$ID %in% rownames(mammalstrait),]
        taxaInfo_mammals<-taxaInfo_mammals[taxaInfo_mammals$ID %in% mammalsID$ID,]        
       
    
        #Birds
        load(file=file.path(data_dir,"birds","birdsID.RData"))
        load(file=file.path(results_dir,"birds","birdstrait.RData"))
        load(file=file.path(results_dir,"birds","50km","occ_birds_list.RData"))
        load(file=file.path(results_dir,"birds","50km","FR_birds.RData"))
        load(file=file.path(data_dir,"birds","taxaInfo_birds.RData"))
        birdsID<-birdsID[birdsID$ID %in% rownames(birdstrait),]
        taxaInfo_birds<-taxaInfo_birds[taxaInfo_birds$ID %in% birdsID$ID,]      

    
#----

        
# Load Phylogeny
    load(file=file.path(data_dir,"mammals","mammalsPhy.RData")) # Both corresponding to Tree1 of the phylogeny
    birdsPhy<-read.tree(file=file.path(data_dir,"birds","birdsPhy.tre")) 
    
# Dropping names not in  ID
    set_mammals <- ape::drop.tip(mammalsPhy,mammalsPhy$tip.label[!is.element(mammalsPhy$tip.label,as.character(gsub(" ", "_", mammalsID$Name)))])
    set_birds <- ape::drop.tip(birdsPhy,birdsPhy$tip.label[!is.element(birdsPhy$tip.label,as.character(gsub(" ", "_", birdsID$Name)))])
    
# Compute Evolutionary Distinctiveness 
    #We use the identical framework (not based on tree but on distance)
    distPhyl_mammals <- cophenetic(set_mammals)
    distPhyl_mammals <- distPhyl_mammals/max(distPhyl_mammals)
    distPhyl_birds <- cophenetic(set_birds)
    distPhyl_birds <- distPhyl_birds/max(distPhyl_birds)
    
    #Square root transformation (following Letten & Cornwell, 2014)
      #allow comparision with function distinctiveness
    distPhyl_mammals <- sqrt(distPhyl_mammals)
    distPhyl_birds <- sqrt(distPhyl_birds)
      
    #Compute ED
    
            # Mammals
                #Evolutionary distinctiveness
                Sim_commu <- matrix(1,1,ncol(distPhyl_mammals))
                colnames(Sim_commu) <- colnames(distPhyl_mammals)
                EDi_mammals <- t(distinctiveness(Sim_commu,distPhyl_mammals))
                colnames(EDi_mammals)<-"EDi"
                
                mammalsID2 = mammalsID
                mammalsID2$Name<-  as.character(gsub(" ", "_", mammalsID2$Name))
                EDi_mammals <- merge(EDi_mammals,mammalsID2, by.x="row.names", by.y="Name" )
                EDi_mammals$EDin <-(EDi_mammals$EDi-min(EDi_mammals$EDi)) / max(EDi_mammals$EDi-min(EDi_mammals$EDi))
                
                #Evolutionary originarluity
                ty <- evol.distinct(set_mammals,type="equal.splits")
                ty <- merge(ty,mammalsID2, by.x="Species", by.y="Name" )
                ty$tyin<-(ty$w-min(ty$w)) / max(ty$w-min(ty$w))
                rownames(ty)<-ty$ID
                ty<- ty[,c(2,4)]
                  
            #Evolutionary distinctiveness
            EDi_mammals$EDin <-(EDi_mammals$EDi-min(EDi_mammals$EDi)) / max(EDi_mammals$EDi-min(EDi_mammals$EDi))
            rownames(EDi_mammals) <- EDi_mammals$ID
            EDi_mammals <- EDi_mammals[,c(2,4)]
            rownames(EDi_mammals) <- EDi_mammals$ID
            
            EDi_mammals <-merge (EDi_mammals,ty, by = "row.names")
            rownames(EDi_mammals) <- EDi_mammals$ID
            EDi_mammals <- EDi_mammals[,c(3,5:7)]
            colnames(EDi_mammals) <- c("EDi", "EDin", "OriDi", "OriDin")
          
            FR_mammals$FR <- merge(FR_mammals$FR,EDi_mammals,by="row.names",all.x = TRUE)
            rownames(FR_mammals$FR) <- FR_mammals$FR$Row.names
            FR_mammals$FR <- FR_mammals$FR[,-1]
            
            # Birds TODO
            Sim_commu <- matrix(1,1,ncol(distPhyl_birds))
            colnames(Sim_commu) <- colnames(distPhyl_birds)
            EDi_birds <- distinctiveness(Sim_commu,distPhyl_birds)
      
            
            data_DR_mammals <- merge(data_DR_mammals,FR_mammals$FR,by ="row.names")
            rownames(data_DR_mammals) <- data_DR_mammals$Row.names
            data_DR_mammals <- data_DR_mammals[,-1]
            
            
            data_DR_mammals <- na.omit(data_DR_mammals)
a <- ggplot(data_DR_mammals, aes(x=DR_class, y=OriDin, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
     geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, 1)) + geom_hline(yintercept=mean(data_DR_mammals$OriDin,na.rm=T),col="red",linetype="dashed") + 
     labs(x = "DR class",y="OriDin")+theme_bw()
            
b <- ggplot(data_DR_mammals, aes(x=DR_class, y=EDin, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, 1)) + geom_hline(yintercept=mean(data_DR_mammals$EDin,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="EDin")+theme_bw()           
            
            
            
            
draw.phylo <- function(FR_data,taxaInfo,set_phylo,taxa) {

  #  FR_data<-FR_mammals
  #  set_phylo <- set_mammals
  #  taxaInfo<- taxaInfo_mammals
  #  taxa="mammals"
  
  
   FR_data<-FR_birds
    set_phylo <- set_birds
    taxaInfo<- taxaInfo_birds
    taxa="birds"
  
      # Create Class DR 
      data_DR<-FR_data$FR
      data_DR$DR_class="NA"
      QD75 <- FR_data$Q$Q75_D
      QD25 <- FR_data$Q$Q25_D
      QR75 <- FR_data$Q$Q75_R
      QR25 <- FR_data$Q$Q25_R
      
      data_DR$DR_class[(data_DR$Din<QD25) & (data_DR$Rin<QR25)]="D25R25"
      data_DR$DR_class[(data_DR$Din>QD75) & (data_DR$Rin>QR75)]="D75R75"
      data_DR$DR_class[(data_DR$Din<QD25) & (data_DR$Rin>QR75)]="D25R75"
      data_DR$DR_class[(data_DR$Din>QD75) & (data_DR$Rin<QR25)]="D75R25"
      data_DR$DR_class[(((data_DR$Din>QD25) & (data_DR$Din<QD75)) & ((data_DR$Rin>QR25) & (data_DR$Rin<QR75)))]="AVG"
      
      data_DR$InvRin=1-data_DR$Rin

      data_DR<-merge(data_DR,taxaInfo,by.x="row.names", by.y="ID")   # ATTENTION
      #data_DR <- data_DR[,-c(11:13,15,16)]
      data_DR <- data.frame(data_DR[,c(1:10)],data_DR$order)
      colnames(data_DR)[11] <- "order"

      
        rownames(data_DR)<-as.character(gsub(" ", "_", data_DR$Name))
        data_DR<-data_DR[rownames(data_DR) %in% set_phylo$tip.label,]
        #order in the same order of the phylo
        data_DR<-data_DR[set_phylo$tip.label,] 
  
        data_DR$cols <- NA
        data_DR$cols[data_DR$DR_class=="D75R75"] <- "red"
    
        
        #The species  Coracina_melas  will cause problems - searching for this species in all the species present will return multiple species, consider renaming! 
        data_DR <- data_DR[which(rownames(data_DR) != "Coracina_melas"),]  
        data_DR <- data_DR[which(rownames(data_DR) != "Lanius_collurio"),] 
        data_DR <- data_DR[which(rownames(data_DR) != "Lanius_excubitor"),]
        data_DR <- data_DR[which(rownames(data_DR) != "Buteo_augur"),]
        data_DR <- data_DR[which(rownames(data_DR) != "Haliaeetus_vocifer"),]
        
        # Prepare order labels
        labelsArc <- na.omit(as.character(unique(data_DR$order)))
        # Change capita in the names of order
        labelsArc %<>% tolower
        labelsArc<-firstup(labelsArc)

        if (taxa=="mammals"){
        labelsArc <- labelsArc[-which(labelsArc == "Dermanoptera")]
        labelsArc <- labelsArc[-which(labelsArc == "Tubulidentata")] 
        labelsArc <- labelsArc[-which(labelsArc == "Microbiotheria") ]} 

        if (taxa=="birds"){
        labelsArc <- labelsArc[-which(labelsArc == "Struthioniformes")]
        labelsArc <- labelsArc[-which(labelsArc == "Caryophyllales")]
        labelsArc <- labelsArc[-which(labelsArc == "Hemiptera")]
        labelsArc <- labelsArc[-which(labelsArc == "Leptosomiformes")]
        labelsArc <- labelsArc[-which(labelsArc == "Opisthocomiformes")]}
        #TODO RESOUDRE CE PROBLEME!# order having only one species and make problem!!!
        
        #Finds nodes of Arc
        nodesArc <- unlist(lapply(labelsArc, function(x){
        node <- phytools::findMRCA(set_phylo,as.character(set_phylo$tip.label[which(as.character(data_DR$order) == x)]), type = "node")
                                   names(node) <- x
                                   node}))
        nodesArc <- nodesArc[order(nodesArc, decreasing = FALSE)]
       
          # Add column binary for Functional Rarity: Yes/no
          rarety <-  data_DR$cols 
          rarety <- as.numeric(as.factor(rarety))
          rarety[is.na(rarety)]<-0
          names(rarety)<-rownames(data_DR)
          
          set_phylo <- ape::drop.tip(set_phylo,set_phylo$tip.label[!is.element(set_phylo$tip.label,as.character(gsub(" ", "_", rownames(data_DR))))])
          
          # plotting PHYLOGENY TREE
          color.terminal.branches(set_phylo, rarety, breaks=4, cols=c("#A6A6A666","red"), edge.width=0.4, show.tip.label=TRUE, non.terminal.col= "#A6A6A666")
          tiplabels(pch = 18, col = data_DR$cols, cex = 0.4 ,offset=5)

          # plotting family labels/arcs
          offset <- rep(c(1.10,1.18,1.26),length(nodesArc)/2)
          for(i in 1:length(nodesArc)){
            
            if(i %in% c(seq(1,length(nodesArc), by = 2))) laboffset <- 0.03
            if(i %in% c(seq(2,length(nodesArc), by = 2))) laboffset <- 0.03
            
            if(i %in% c(1,4,7,10,13,16,19,22,25,28,31,34,37)){  #If odd 
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
              
            if(i %in% c(2,5,8,11,14,17,20,23,26,29,32,35)){  #If odd 
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
                
            if(i %in% c(3,6,9,12,15,18,21,24,27,30,33,36)){  #If odd 
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
          if (taxa=="mammals") text(x=270,y=130-(i*10),labels=paste0(i,": ",names(nodesArc)[i]),cex=0.4)
          if (taxa=="birds"){text(x=170,y=130-(i*5),labels=paste0(i,": ",names(nodesArc)[i]),cex=0.4) 
           } 
           }
          
#---
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
                  
                      #collapse or resolve multichotomies in phylogenetic trees TODO check that is mean exactely because need it
                      set_tree<-di2multi(set_tree)
                      
                      #Compute D and statistic
                      FR_PhyloD <- comparative.data(set_tree, data_DR,"species",na.omit=FALSE)
                      FR_PhyloD <- phylo.d(FR_PhyloD, binvar=rarety,permut=permut)
                      
                      #The estimated D value
                      estimated_D <- FR_PhyloD$DEstimate
                      #A p value,giving the result of testing whether D is significantly different from one
                      Pval1 <- FR_PhyloD$Pval1
                      #A p value, giving the result of testing whether D is significantly different from zero
                      Pval0 <- FR_PhyloD$Pval0
                      
                      Dstat <- data.frame(estimated_D,Pval1,Pval0)
                      return(Dstat)
                      },mc.cores= proc)
  
          }
#Mammals
D_mammals <- do.call(rbind,D.phylogeny(ids=1:100,proc=3,data_DR=data_DR,taxa="mammals",permut=1000))
save(D_mammals,file=file.path(results_dir,"mammals","50km","D_mammals.RData"))
D_mammals_plot<-ggplot(D_mammals, aes(estimated_D)) + geom_density(adjust = 1.5,alpha = 0.1,fill="red",colour="red") + xlim(0, 1)+theme_bw()+  labs(x = "D")+
  theme(axis.title=element_text(size=8),axis.text.x = element_text(size=6))
D_mammals_plot<-print(D_mammals_plot, vp=viewport(.5, .5, .17, .15))

#birds
D_birds <- do.call(rbind,D.phylogeny(ids=1:100,proc=3,data_DR=data_DR,taxa="birds",permut=1000))
save(D_birds,file=file.path(results_dir,"birds","50km","D_birds.RData"))
D_birds_plot<-ggplot(D_birds, aes(estimated_D)) + geom_density(adjust = 1.5,alpha = 0.1,fill="red",colour="red") + xlim(0, 1)+theme_bw()+  labs(x = "D")+
  theme(axis.title=element_text(size=8))
D_birds_plot<-print(D_birds_plot, vp=viewport(.5, .5, .17, .15))# vp=viewport(.12, .85, .24, .22))


