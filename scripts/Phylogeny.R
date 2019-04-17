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
library(viridis)
library(caper)

#LOAD TRAITS MAPS AND DISTRIB----

    # Load traits and distrib 

              #Mammals
              load(file=file.path(results_dir,"mammals","mammalsID.RData"))
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
  load(file=file.path(data_dir,"mammals","mammalsPhy.rdata")) # Both corresponding to Tree1 of the phylogeny
  birdsPhy<-read.tree(file=file.path(data_dir,"birds","birdsPhy.tre")) 
  i <- 3
  birdsPhy<-read.tree(file=file.path(data_dir,"birds","alltrees", paste0("BirdzillaHackett1_",i,".tre")))

          # Dropping names not in  ID
            set_mammals <- ape::drop.tip(mammalsPhy,mammalsPhy$tip.label[!is.element(mammalsPhy$tip.label,as.character(gsub(" ", "_", mammalsID$Name)))])
            set_birds <- ape::drop.tip(birdsPhy,birdsPhy$tip.label[!is.element(birdsPhy$tip.label,as.character(gsub(" ", "_", birdsID$Name)))])


#SOUS JEUX DE DONNEE POUR FAIRE DES TEST 
#set_birds <- ape::drop.tip(birdsPhy,birdsPhy$tip.label[!is.element(birdsPhy$tip.label,sample(as.character(gsub(" ", "_", birdsID$Name)),200))])

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
          #Prepare data: DR class for each species
              #Mammals
              load(file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))
              data_DR_mammals<-FR_mammals$FR
              QD75 <- FR_mammals$Q$Q75_D
              QD25 <- FR_mammals$Q$Q25_D
              QR75 <- FR_mammals$Q$Q75_R
              QR25 <- FR_mammals$Q$Q25_R
              
              data_DR_mammals$DR_class[(data_DR_mammals$Din<QD25) & (data_DR_mammals$Rin<QR25)]="D25R25"
              data_DR_mammals$DR_class[(data_DR_mammals$Din>QD75) & (data_DR_mammals$Rin>QR75)]="D75R75"
              data_DR_mammals$DR_class[(data_DR_mammals$Din<QD25) & (data_DR_mammals$Rin>QR75)]="D25R75"
              data_DR_mammals$DR_class[(data_DR_mammals$Din>QD75) & (data_DR_mammals$Rin<QR25)]="D75R25"
              data_DR_mammals$DR_class[(((data_DR_mammals$Din>QD25) & (data_DR_mammals$Din<QD75)) & ((data_DR_mammals$Rin>QR25) & (data_DR_mammals$Rin<QR75)))]="AVG"
              
              data_DR_mammals<-data.frame(data_DR_mammals[,"DR_class"],row.names = rownames(data_DR_mammals))
              colnames(data_DR_mammals) <- "DR_class"
              
              #Birds
              data_DR_birds<-FR_birds$FR
              QD75 <- FR_birds$Q$Q75_D
              QD25 <- FR_birds$Q$Q25_D
              QR75 <- FR_birds$Q$Q75_R
              QR25 <- FR_birds$Q$Q25_R
              
              data_DR_birds$DR_class[(data_DR_birds$Din<QD25) & (data_DR_birds$Rin<QR25)]="D25R25"
              data_DR_birds$DR_class[(data_DR_birds$Din>QD75) & (data_DR_birds$Rin>QR75)]="D75R75"
              data_DR_birds$DR_class[(data_DR_birds$Din<QD25) & (data_DR_birds$Rin>QR75)]="D25R75"
              data_DR_birds$DR_class[(data_DR_birds$Din>QD75) & (data_DR_birds$Rin<QR25)]="D75R25"
              data_DR_birds$DR_class[(((data_DR_birds$Din>QD25) & (data_DR_birds$Din<QD75)) & ((data_DR_birds$Rin>QR25) & (data_DR_birds$Rin<QR75)))]="AVG"
              
              data_DR_birds<-data.frame(data_DR_birds[,"DR_class"],row.names = rownames(data_DR_birds))
              colnames(data_DR_birds) <- "DR_class"  
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
          EDi_mammals$EDin <-(EDi_mammals$EDi-min(EDi_mammals$EDi)) / max(EDi_mammals$EDi-min(EDi_mammals$EDi))
          rownames(EDi_mammals) <- EDi_mammals$ID
          
          #Evolutionary originarlity
          ty <- evol.distinct(set_mammals,type="equal.splits")
          ty <- merge(ty,mammalsID2, by.x="Species", by.y="Name" )
          ty$tyin<-(ty$w-min(ty$w)) / max(ty$w-min(ty$w))
          rownames(ty)<-ty$ID
          ty<- ty[,c(2,4)]
          
          #Merge Data
          EDi_mammals <-merge (EDi_mammals,ty, by = "row.names")
          rownames(EDi_mammals) <- EDi_mammals$ID
          EDi_mammals <- EDi_mammals[,c(3,5:7)]
          colnames(EDi_mammals) <- c("EDi", "EDin", "OriDi", "OriDin")
          
          FR_mammals$FR <- merge(FR_mammals$FR,EDi_mammals,by="row.names",all.x = TRUE)
          rownames(FR_mammals$FR) <- FR_mammals$FR$Row.names
          FR_mammals$FR <- FR_mammals$FR[,-1]
          
          data_DR_mammals <- merge(data_DR_mammals,FR_mammals$FR,by ="row.names")
          rownames(data_DR_mammals) <- data_DR_mammals$Row.names
          data_DR_mammals <- data_DR_mammals[,-1]


  # Birds
          
          #Evolutionary distinctiveness
          Sim_commu <- matrix(1,1,ncol(distPhyl_birds))
          colnames(Sim_commu) <- colnames(distPhyl_birds)
          EDi_birds <- t(distinctiveness(Sim_commu,distPhyl_birds))
          colnames(EDi_birds)<-"EDi"
          
          birdsID2 = birdsID
          birdsID2$Name<-  as.character(gsub(" ", "_", birdsID2$Name))
          EDi_birds <- merge(EDi_birds,birdsID2, by.x="row.names", by.y="Name" )
          EDi_birds$EDin <-(EDi_birds$EDi-min(EDi_birds$EDi)) / max(EDi_birds$EDi-min(EDi_birds$EDi))
          EDi_birds$EDin <-(EDi_birds$EDi-min(EDi_birds$EDi)) / max(EDi_birds$EDi-min(EDi_birds$EDi))
          rownames(EDi_birds) <- EDi_birds$ID
          
          #Evolutionary originarlity
          ty <- evol.distinct(set_birds,type="equal.splits")
          ty <- merge(ty,birdsID2, by.x="Species", by.y="Name" )
          ty$tyin<-(ty$w-min(ty$w)) / max(ty$w-min(ty$w))
          rownames(ty)<-ty$ID
          ty<- ty[,c(2,4)]
          
          #Merge Data
          EDi_birds <-merge (EDi_birds,ty, by = "row.names")
          rownames(EDi_birds) <- EDi_birds[,1]
          EDi_birds <- EDi_birds[,c(3,5:7)]
          colnames(EDi_birds) <- c("EDi", "EDin", "OriDi", "OriDin")
          
          FR_birds$FR <- merge(FR_birds$FR,EDi_birds,by="row.names",all.x = TRUE)
          rownames(FR_birds$FR) <- FR_birds$FR$Row.names
          FR_birds$FR <- FR_birds$FR[,-1]
          
          data_DR_birds <- merge(data_DR_birds,FR_birds$FR,by ="row.names")
          rownames(data_DR_birds) <- data_DR_birds$Row.names
          data_DR_birds <- data_DR_birds[,-1]

#Plot Evolutionary originality and distinctiveness  
               
          
          data_DR_mammalsPlot <- na.omit(data_DR_mammals)
          data_DR_birdsPlot <- na.omit(data_DR_birds)
          a <- ggplot(data_DR_mammalsPlot, aes(x=DR_class, y=OriDin, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + #scale_fill_manual(values=col_br)+
            geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, 1)) + geom_hline(yintercept=mean(data_DR_mammalsPlot$OriDin,na.rm=T),col="red",linetype="dashed") + 
            labs(x = "DR class",y="Evolutionary originarlity")+theme_bw()
          
          b <- ggplot(data_DR_mammalsPlot, aes(x=DR_class, y=EDin, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + #scale_fill_manual(values=col_br)+
            geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, 1)) + geom_hline(yintercept=mean(data_DR_mammalsPlot$EDin,na.rm=T),col="red",linetype="dashed") + 
            labs(x = "DR class",y="Evolutionary distinctiveness")+theme_bw()           
          
          
          c <- ggplot(data_DR_birdsPlot, aes(x=DR_class, y=OriDin, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + #scale_fill_manual(values=col_br)+
            geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, 1)) + geom_hline(yintercept=mean(data_DR_birdsPlot$OriDin,na.rm=T),col="red",linetype="dashed") + 
            labs(x = "DR class",y="Evolutionary originarlity")+theme_bw()
            
          
          d <- ggplot(data_DR_birdsPlot, aes(x=DR_class, y=EDin, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + #scale_fill_manual(values=col_br)+
            geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, 1)) + geom_hline(yintercept=mean(data_DR_birdsPlot$EDin,na.rm=T),col="red",linetype="dashed") + 
            labs(x = "DR class",y="Evolutionary distinctiveness")+theme_bw()   
          #Normal que le EDin fonctionne car m??me construction donc ?? ne pas utilis??  

          gridExtra::grid.arrange(a,b,c,d,ncol= 2)  
draw.phylo <- function(FR_data,taxaInfo,set_phylo,taxa) {
  
  FR_data<-FR_mammals
  set_phylo <- set_mammals
  taxaInfo<- taxaInfo_mammals
  taxa="mammals"
  
  #FR_data<-FR_birds
  #set_phylo <- set_birds
  #taxaInfo<- taxaInfo_birds
  #taxa="birds"
  
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
  
  rownames(data_DR)<-as.character(gsub(" ", "_", data_DR$Name))
  data_DR<-data_DR[rownames(data_DR) %in% set_phylo$tip.label,]
  #order in the same order of the phylo
  data_DR<-data_DR[set_phylo$tip.label,] 
  
  data_DR$colsD75R75 <- NA
  data_DR$colsD75R75[data_DR$DR_class=="D75R75"] <- "orangered"
  
  data_DR$colsD25R25 <- NA
  data_DR$colsD25R25[data_DR$DR_class=="D25R25"] <- "#E7B800"
  
  data_DR$colsAVG <- NA
  data_DR$colsAVG[data_DR$DR_class=="AVG"] <- "#00AFBB"
  
  
  #The species  Coracina_melas  will cause problems - searching for this species in all the species present will return multiple species, consider renaming! 
  data_DR <- data_DR[which(rownames(data_DR) != "Coracina_melas"),]  
  data_DR <- data_DR[which(rownames(data_DR) != "Lanius_collurio"),] 
  data_DR <- data_DR[which(rownames(data_DR) != "Lanius_excubitor"),]
  data_DR <- data_DR[which(rownames(data_DR) != "Buteo_augur"),]
  data_DR <- data_DR[which(rownames(data_DR) != "Haliaeetus_vocifer"),]
  
  # Change capita in the names of order
  data_DR$order %<>% tolower
  data_DR$order<-firstup(data_DR$order)
  
  # Prepare order labels
  labelsArc <- na.omit(as.character(unique(data_DR$order)))
  
  if (taxa=="mammals"){
    labelsArc <- labelsArc[-which(labelsArc == "Dermoptera")]
    labelsArc <- labelsArc[-which(labelsArc == "Tubulidentata")] 
    labelsArc <- labelsArc[-which(labelsArc == "Microbiotheria") ]} 
  
  if (taxa=="birds"){
    labelsArc <- labelsArc[-which(labelsArc == "Struthioniformes")]
    labelsArc <- labelsArc[-which(labelsArc == "Caryophyllales")]
    labelsArc <- labelsArc[-which(labelsArc == "Hemiptera")]
    labelsArc <- labelsArc[-which(labelsArc == "Leptosomiformes")]
    labelsArc <- labelsArc[-which(labelsArc == "Opisthocomiformes")]}
  #TODO RESOUDRE CE PROBLEME!# order having only one species and make problem!!!
  #A VIRER
  #Finds nodes of Arc
  nodesArc <- unlist(lapply(labelsArc, function(x){
    node <- phytools::findMRCA(set_phylo,as.character(set_phylo$tip.label[which(as.character(data_DR$order) == x)]), type = "node")
    names(node) <- x
    node}))
 

  nodesArc <- nodesArc[order(nodesArc, decreasing = FALSE)]

  # Add column binary for Functional Rarity: Yes/no
  rarety <-  data_DR$colsD75R75 
  rarety <- as.numeric(as.factor(rarety))
  rarety[is.na(rarety)]<-0
  names(rarety)<-rownames(data_DR)

    
  set_phylo <- ape::drop.tip(set_phylo,set_phylo$tip.label[!is.element(set_phylo$tip.label,as.character(gsub(" ", "_", rownames(data_DR))))])
  data_DR <-data_DR[as.character(gsub(" ", "_", rownames(data_DR)))%in%set_phylo$tip.label,]
  # plotting PHYLOGENY TREE
  color.terminal.branches(set_phylo, rarety, breaks=4, cols=c("#A6A6A63F","orangered"), edge.width=0.4, show.tip.label=TRUE,non.terminal.col= "#A6A6A63F")
  tiplabels(pch = 16, col = data_DR$colsD75R75, cex = 0.4 ,offset=9)
  tiplabels(pch = 16, col = data_DR$colsD25R25, cex = 0.4 ,offset=6)
  tiplabels(pch = 16, col = data_DR$colsAVG, cex = 0.4 ,offset=2)
  
  
  # plotting family labels/arcs
  offset <- rep(c(1.15,1.23),length(nodesArc))
  
  
  for(i in 1:length(nodesArc)){
    
    if(i %in% c(seq(1,length(nodesArc), by = 2))) laboffset <- 0.03
    if(i %in% c(seq(2,length(nodesArc), by = 2))) laboffset <- 0.03
    
   if(i %in% c(1,4,7,10,13,16,19,22,25,28,31,34,37)){  #If odd 
      arc.cladelabels(text= paste0(i,""), #paste0(names(nodesArc)[i])
                      node=nodesArc[i],
                      ln.offset=offset[i],
                      lab.offset=offset[i]+laboffset, 
                      cex = 1, 
                      colarc = "gray82",
                      lwd = 1, 
                      lty = 1, cextext=0.5,
                      orientation = "curved", mark.node = FALSE,
                      col="gray82",coltext="gray82")}
    
    if(i %in% c(2,5,8,11,14,17,20,23,26,29,32,35)){  #If odd 
      arc.cladelabels(text= paste0(i,""), #paste0(names(nodesArc)[i])
                      node=nodesArc[i],
                      ln.offset=offset[i],
                      lab.offset=offset[i]+laboffset, 
                      cex = 1, 
                      colarc = "gray75",
                      lwd = 1, 
                      lty = 1, cextext=0.5,
                      orientation = "curved",
                      mark.node = FALSE,col="gray75",coltext="gray75")}
    
    if(i %in% c(3,6,9,12,15,18,21,24,27,30,33,36)){  #If odd 
      arc.cladelabels(text= paste0(i,""), #paste0(names(nodesArc)[i])
                      node=nodesArc[i],
                      ln.offset=offset[i],
                      lab.offset=offset[i]+laboffset, 
                      cex = 1, 
                      colarc = "gray68",
                      lwd = 1, 
                      lty = 1, cextext=0.5,
                      orientation = "curved",
                      mark.node = FALSE,col="gray68",coltext="gray68")}
    
    
    #Add the names of Order
    if (taxa=="mammals") text(x=270,y=130-(i*10),labels=paste0(i,": ",names(nodesArc)[i]),cex=0.5)
    if (taxa=="birds")   text(x=170,y=130-(i*5),labels=paste0(i,": ",names(nodesArc)[i]),cex=0.4) 
  }
}

#---
#Compute FRITZ to know if functional rare species are packaged      

data_DR$rarety <- data_DR$colsD75R75
data_DR$rarety <- as.numeric(as.factor(data_DR$rarety))
data_DR$rarety[is.na(data_DR$rarety)]<-0



data_DR$D75R75 <- data_DR$colsD75R75
data_DR$D75R75 <- as.numeric(as.factor(data_DR$D75R75))
data_DR$D75R75[is.na(data_DR$D75R75)]<-0

data_DR$D25R25 <- data_DR$colsD25R25
data_DR$D25R25 <- as.numeric(as.factor(data_DR$D25R25))
data_DR$D25R25[is.na(data_DR$D25R25)]<-0

data_DR$AVG <- data_DR$colsAVG
data_DR$AVG <- as.numeric(as.factor(data_DR$AVG))
data_DR$AVG[is.na(data_DR$AVG)]<-0

#Should be done with the 100 trees of mammals et 100 trees of birds
data_DR$species<-rownames(data_DR)

D.phylogeny <- function(ids,proc,data_DR,taxa,permut) {
  #proc <- 2
  #taxa <- "birds"
  #data_DR <- data_DR
  #permut <- 10
  #ids <- 1:2
  #indice="D75R75"
  
  mclapply(ids,function(id) { 
    if (taxa == "mammals") tree<-read.tree(file=file.path(data_dir,"mammals","alltrees", paste0("FritzTree_mammals_updateCarnivora2012DEF",id,".tre"))) 
    
    if (taxa == "birds")   tree<-read.tree(file=file.path(data_dir,"birds","alltrees", paste0("BirdzillaHackett1_",id,".tre"))) 
    
    set_tree<- ape::drop.tip(tree,tree$tip.label[!is.element(tree$tip.label,as.character(gsub(" ", "_", rownames(data_DR))))])
    set_tree$node.label <- NULL
    
    #collapse or resolve multichotomies in phylogenetic trees TODO check that is mean exactely because need it
    set_tree<-di2multi(set_tree)
    
    #Compute D and statistic
    FR_PhyloD <- caper::comparative.data(set_tree, data_DR,"species",na.omit=FALSE)
    FR_PhyloD <- caper::phylo.d(FR_PhyloD, binvar=D75R75,permut=permut)
    
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
#D_mammalsAVG <- do.call(rbind,D.phylogeny(ids=1:100,proc=25,data_DR=data_DR,taxa="mammals",permut=1000,indice="D25R25"))
#save(D_mammalsAVG,file=file.path(results_dir,"mammals","50km","D_mammalsAVG.RData"))
#D_mammals <- do.call(rbind,D.phylogeny(ids=1:100,proc=3,data_DR=data_DR,taxa="mammals",permut=1000))
#save(D_mammals,file=file.path(results_dir,"mammals","50km","D_mammals.RData"))
#D_mammalsD25R25 <- do.call(rbind,D.phylogeny(ids=1:100,proc=25,data_DR=data_DR,taxa="mammals",permut=1000))
#save(D_mammalsD25R25,file=file.path(results_dir,"mammals","50km","D_mammalsD25R25.RData"))
D_mammalsD75R75 <- do.call(rbind,D.phylogeny(ids=1:100,proc=25,data_DR=data_DR,taxa="mammals",permut=1000))
save(D_mammalsD75R75 ,file=file.path(results_dir,"mammals","50km","D_mammalsD75R75.RData"))

load(file=file.path(results_dir,"mammals","50km","D_mammalsD75R75.RData"))
load(file=file.path(results_dir,"mammals","50km","D_mammalsAVG.RData"))
load(file=file.path(results_dir,"mammals","50km","D_mammalsD25R25.RData"))

D_all <- data.frame(rbind(D_mammalsD75R75,D_mammalsAVG,D_mammalsD25R25),DR_class=c(rep("D75R75",100),rep("AVG",100),rep("D25R25",100)))
D_all_mammals_plot<-ggplot(D_all, aes(estimated_D,color=DR_class,fill=DR_class)) + geom_density(adjust = 1.5,alpha = 0.1) + xlim(0, 1)+theme_bw()+  labs(x = "D")+
  scale_color_manual(values = c("#00AFBB", "#E7B800","orangered")) + scale_fill_manual(values = c("#00AFBB", "#E7B800","orangered")) +
  theme(axis.title=element_text(size=8),axis.text.x = element_text(size=6))+theme_bw()
  D_all_mammals_plot<-print(D_all_mammals_plot, vp=viewport(.5, .5, .17, .15))
  



z
#birds
D_birdsD75R75 <- do.call(rbind,D.phylogeny(ids=1:100,proc=25,data_DR=data_DR,taxa="birds",permut=1000))
save(D_birdsD75R75,file=file.path(results_dir,"birds","50km","D_birds.RData"))
#D_birdsAVG <- do.call(rbind,D.phylogeny(ids=1:100,proc=25,data_DR=data_DR,taxa="birds",permut=1000))
#save(D_birdsAVG,file=file.path(results_dir,"birds","50km","D_birdsAVG.RData"))
#D_birdsD25R25 <- do.call(rbind,D.phylogeny(ids=1:100,proc=25,data_DR=data_DR,taxa="birds",permut=1000))
#save(D_birdsD25R25,file=file.path(results_dir,"birds","50km","D_birdsD25R25.RData"))


load(file=file.path(results_dir,"birds","50km","D_birds.RData"))
load(file=file.path(results_dir,"birds","50km","D_birdsAVG.RData"))
load(file=file.path(results_dir,"birds","50km","D_birdsD25R25.RData"))
D_birds_plot<-ggplot(D_birds, aes(estimated_D)) + geom_density(adjust = 1.5,alpha = 0.1,fill="red",colour="red") + xlim(0, 1)+theme_bw()+  labs(x = "D")+
  theme(axis.title=element_text(size=8))
D_birds_plot<-print(D_birds_plot, vp=viewport(.5, .5, .17, .15))# vp=viewport(.12, .85, .24, .22))

D_birds_plot<-ggplot(D_birds, aes(estimated_D)) + geom_histogram(alpha = 0.1,fill="red",colour="red") + xlim(0.25, 0.75)+theme_bw()+  labs(x = "D")+
  theme(axis.title=element_text(size=8))

