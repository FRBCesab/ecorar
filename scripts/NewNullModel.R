
library(funrar)
library(moments)
library(ggplot2)
library(parallel)
library(ade4)
library(dplyr)
library(gridExtra)
library(cluster)
library(rgdal)
library(plyr)



load(file=file.path(results_dir,"mammals/mammalsID.RData"))
load(file=file.path(results_dir,"mammals","mammalstrait.RData"))
load(file=file.path(results_dir,"mammals","50km","occ_mammals_list.RData"))
mammalsID<-mammalsID[mammalsID$ID %in% rownames(mammalstrait),]

load(file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))
load(file=file.path(results_dir,"birds","50km","FR_birds.RData"))




Ri<-data.frame(table(unlist(occ_mammals_list))/length(occ_mammals_list))
rownames(Ri)<-Ri[,1]
Ri<-Ri[,-1, drop = FALSE]
Ri<-1-Ri
colnames(Ri)<-"Rin"
Rin<-(Ri-min(Ri)) / max(Ri-min(Ri))



Q75_D <- FR_mammals$Q$Q75_D
Q75_R <- FR_mammals$Q$Q75_R



# Step 1 : Create matrice for each species 1000 column of class


#Randomize names of species in functional trait matrix

mammalstrait_random <- mammalstrait
rownames(mammalstrait_random) <-  sample (row.names (mammalstrait))

          ###Format the traits 
          diet <- prep.fuzzy(mammalstrait_random[,1:10], col.blocks = ncol(mammalstrait_random[,1:10]), label = "diet")
          
          ForStrat <- mammalstrait_random$ForStrat %>% as.data.frame()
          row.names(ForStrat) <- row.names(mammalstrait_random)
          
          Activity <- apply(mammalstrait_random[,12:14],2,as.numeric) %>% as.data.frame()
          row.names(Activity) <- row.names(mammalstrait_random)
          Activity=prep.binary(Activity, col.blocks = ncol(Activity), label = "Activity")
          
          bodymass <- log(mammalstrait_random$BodyMass.Value) %>% as.data.frame()
          row.names(bodymass) <- row.names(mammalstrait_random)
          bodymass <- bodymass/(max(bodymass, na.rm=T)-min(bodymass, na.rm=T))
          
          ###Compute the dist matrix
          disTraits_mammals <- dist.ktab(ktab.list.df(list(diet, ForStrat, Activity, bodymass)), c("F","N","B","Q"), scan = FALSE) %>% as.matrix()

          ####Compute Di (Global distinctiveness)
          Sim_commu <- matrix(1,1,ncol(disTraits_mammals))
          colnames(Sim_commu) <- colnames(disTraits_mammals)
          Din_random<-t(distinctiveness(Sim_commu,disTraits_mammals))
          colnames(Din_random)<-"Din"
          Din_random<-(Din_random-min(Din_random)) / max(Din_random-min(Din_random))
          
#Random names of occurence
          Rin_random <-  Rin
          #rownames(Rin_random) <-  sample(row.names(Rin))
          
          FR_data_random<-data.frame(Din_random,Rin_random)

          #Associate with group       
          subD75R75 <- mclapply(ids,function(id) {
            spe_sub <- occ_mat_list[[id]]
            spe_sub_D <- spe_sub[FR_data_random[spe_sub,"Din"]>Q75_D]
            spe_sub_D[FR_data_random[spe_sub_D,"Rin"]>Q75_R]
            
          },mc.cores = proc)
          names(subD75R75) <- ids   
          
          data.frame(unlist(subD75R75))
          