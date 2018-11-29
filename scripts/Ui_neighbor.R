rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=TRUE,who="NL")

#----
#Library
library(gridExtra)
library(funrar)
library(moments)
library(ggplot2)
library(parallel)
library(ade4)
library(dplyr)
library(gridExtra)
library(cluster)

#----
#Load DATA
    # Load spatial grid for plotting
    map<-readOGR(file.path(data_dir,"ReferenceGrid10Km","gridLand10km.B.shp"))
    #names of each cell
    ID_cell<-rownames(map@data)
    # Matrice of neighbor
    load(file.path(data_dir,"ReferenceGrid10Km","Mat_neighbour.RData"))# names of the file  = names(occ)-1
    names(Mat_neighbour)<-rownames(occ_mammals)
    
    # Load traits and distrib 
    load(file=file.path(results_dir,"mammals/mammalsID.RData"))
    load(file=file.path(results_dir,"mammals/mammalstrait.RData"))
    load(file=file.path(data_dir,"mammals/occ_mammals_sparseM.RData"))
    
    #Commun ID for mammalsID/occ_mammals/traitmammals ---
    mammalsID<-mammalsID[mammalsID$checkname %in% mammalstrait$checkname,]
    mammalstrait<-merge(mammalsID,mammalstrait,by="checkname")
    rownames(mammalstrait)<-mammalstrait$ID
    mammalstrait<-mammalstrait[,-c(2,3)]
    occ_mammals <- occ_mammals[,colnames(occ_mammals)  %in% mammalsID$ID]
    
    #Select species with at least one occurence
    occ_mammals <- occ_mammals[,colSums(occ_mammals)>0]
    
    load(file=file.path(results_dir,"mammals/disTraits_mammals.RData"))

    load(file=file.path(data_dir,"mammals/occ_mammals_list.RData"))
    occ_mammals_list<-MammalPresence
    rm(MammalPresence)
    
    #Work with a list, delete species that are not in mammalsID
    for (i in 1: length(occ_mammals_list)){
      occ_mammals_list[[i]]<- occ_mammals_list[[i]][occ_mammals_list[[i]]%in%mammalsID$ID]
      print(paste0('i',i))
    }
    save(occ_mammals_list,file=file.path(results_dir,"mammals/occ_mammals_list.RData"))

#----
#Function to compute Ui of each species inside cell where it is occuring, allow consideration of neighbour cells
Ui.funk<-function(occ_mat_list,sp,dist_traits,mat_neigh,proc) {                
  
  # occ_mat_list=occ_mammals_list
  # sp=  "sp140" 
  # dist_traits=disTraits_mammals
  # mat_neigh=Mat_neighbour
  # proc=4
  
  #Select cell with cell where species "sp" is present
  
  cell<-lapply(lapply(occ_mat_list, function(x){x==sp}),sum)>0
  cell<-cell[cell=="TRUE"]

  
  occ_mat_list[[lapply(lapply(occ_mat_list, function(x){x==sp}),sum)>0]]
  
  
  names(occ_mat_list%in%sp)
  cell<-occ_mat_list[,sp]
  
  
  cell<-cell[cell>0]
  ids_cell<-names(cell)

  #id<-ids_cell
        if (length(ids_cell)<proc) proc_real=1 else proc_real=proc
        
        Ui_cell <- do.call(rbind,mclapply(ids_cell,function(id){ 
          spot<-as.data.frame(occ_mat[id,])
          spot<-t(subset(spot,spot[,1]>0))
          
          if(sum(spot)==1){Ui_Sp_cell<-NA
          
          } else{
            disT_cell<-dist_traits[rownames(dist_traits)%in%colnames(spot),]
            disT_cell<-disT_cell[,colnames(disT_cell)%in%colnames(spot)]
            FR_cell <- funrar(spot, disT_cell, rel_abund = FALSE)
            Ui_Sp_cell<-FR_cell$Ui[FR_cell$Ui$species%in%sp,]$Ui
          } 
          
    #Compute Ui considering cell and neigbor
        spot_neigh<-mat_neigh[[id]]

    #If the cell haven't neighbor
        if(length(spot_neigh)==0){Ui_Sp_neigh<-NA
    
    }else{
      spot_neigh<-as.matrix(occ_mat[rownames(occ_mat)%in%c(id,spot_neigh),])
      spot_neigh<-spot_neigh[,apply(spot_neigh,2,sum)>0]
      #If the species is alone 
      if(dim(data.frame(spot_neigh))[2]==1) {Ui_Sp_neigh<-NA
      }else{
        disT_cell_neigh <- dist_traits[rownames(dist_traits)%in%colnames(spot_neigh),]
        disT_cell_neigh <- disT_cell_neigh[,colnames(disT_cell_neigh)%in%colnames(spot_neigh)]
        
        FR_cell_neigh <- funrar(spot_neigh, disT_cell_neigh, rel_abund = FALSE)
        Ui_Sp_neigh<-FR_cell_neigh$Ui[FR_cell_neigh$Ui$species%in%sp,]$Ui
      }  
    }
    
    res <- cbind.data.frame(Ui_Sp_cell,Ui_Sp_neigh)
    
    names(res) <- c('Ui_Sp_cell','Ui_Sp_neigh')
    
    return(res)},mc.cores = proc_real))
}


#----
#mammals          
varSP <-  colnames(occ_mammals)

####... LOOP
Ui_mammals = list()
ptm <- proc.time()
    for(i in 1:length(varSP))
    {
      print(paste0("i=",i, "/sp=",varSP[i]," / cells= ",sum(occ_mammals[,varSP[i]])," / " ,length(varSP)," mammals"))
      #print(paste0("sp=",i," / cells= ",sum(occ_mammals[,varSP[i]])," / " ,length(varSP)," mammals"))
      LM <- Ui.funk(occ_mat=occ_mammals,sp=varSP[i],dist_traits=disTraits_mammals,mat_neigh=Mat_neighbour,proc=4)
      Ui_mammals[[length(Ui_mammals)+1]] <- LM
      names(Ui_mammals)[[length(Ui_mammals)]] <- varSP[i]
    }
    proc.time()-ptm
save(Ui_mammals, file=file.path(results_dir,"mammals/Ui_mammals.RData"))

Ui_mammals<-do.call(rbind.data.frame, lapply(Ui_mammals, function(x) colMeans(x[1:2])))
save(Ui_mammals, file=file.path(results_dir,"mammals/Ui_mammals.RData"))
