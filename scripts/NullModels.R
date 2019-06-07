#Null model
rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=TRUE,who="NM")

library(parallel)
library(plyr)
library(vegan)
library(qdapTools)


#Prepare data: DR class for each species
load(file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))
#Mammals
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


load(file=file.path(results_dir,"birds","50km","FR_birds.RData"))        
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

##SECOND VERSION OF NULL MODEL (SUGGESTED BY ANNETTE)

#Random distribution but keep the DR_class for each species


# Transform in matrix

load(file=file.path(results_dir,"mammals","50km","occ_mammals_list.RData"))
load(file=file.path(results_dir,"mammals","data_DR_mammals.RData"))
load(file=file.path(results_dir,"mammals","50km","funk_mammals.RData"))


occ_mammals_mat<-mtabulate(occ_mammals_list)!=0
occ_mammals_mat<-occ_mammals_mat[apply(occ_mammals_mat,1,sum)>0,]
#TEst petit jeux de données
#funk_mammals[funk_mammals$D75R75==12,]
#occ_mammals_mat <- occ_mammals_mat[rownames(occ_mammals_mat)%in% as.character(158899:158910),]
#occ_mammals_mat<- occ_mammals_mat[,colnames(occ_mammals_mat) %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D75R75"))]

for (i in 1:ncol(occ_mammals_mat)){
  occ_mammals_mat[,i] <- as.numeric(occ_mammals_mat[,i])
}


rep<-200

simu <- mclapply(1:rep,function(i){
  nullmod <- nullmodel(occ_mammals_mat,method="curveball")
  sim_matrix <- simulate(nullmod, nsim=1)
  colnames(sim_matrix) <-  colnames(occ_mammals_mat)
  sim_matrix <- sim_matrix[,colnames(sim_matrix) %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D75R75")),]
  return(sim_matrix)
},mc.cores=15)


Null_res <-mclapply(1:length(simu),function(i){
  Null_mean <- data.frame(mean=apply(simu[[i]],1,mean),sd=apply(simu[[i]],1,sd))
  },mc.cores=15)
  
Null_mean <- apply(data.frame(lapply(Null_res,function(x) x[,1])),1,mean)
Null_sd <- apply(data.frame(lapply(Null_res,function(x) x[,2])),1,mean)

funk_mammals<-funk_mammals[funk_mammals$TD_sp>0,]
SES_total_mammals <- data.frame(cell=funk_mammals$cell, D75R75 = (funk_mammals$D75R75 - Null_mean)/Null_sd)
save(SES_total_mammals,file=file.path(results_dir,"mammals","50km","SES_total_mammals_Swap.RData"))
boxplot(SES_total_mammals$D75R75)




#THIRD VERSION OF NULL MODEL (SUGGESTED BY BRYAN E. )
#Random the functional trait but keep the rarity to know if FR is link to R only


###FIRST VERSION OF NULL MODEL
#RANDOM DISTRIBUTION OF THE DR CLASS

##Generate number of species per DR_class
Nb.DR_class<- function(ids,proc,occ_mat_list,data_DR_null){
  
  # proc <- 50
  # occ_mat_list <- occ_mammals_list
  # ids <- names(occ_mat_list)
  # data_DR_null <- data_DR_randomize_mammals[[i]]
  
  Null_funk_all <-do.call(rbind,mclapply(ids,function(id) {    #cat("id:",i,"\n")
    
    #id <- ids[1]
    
    #Subset the species present in the site i 
    spe_sub <- unique(occ_mat_list[[id]])
    
    #Number of species within a sites with FR values above quantiles of species and functional distribution 
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="D75R75"))])>0) D75R75=length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="D75R75"))]) else D75R75=0
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="D25R75"))])>0) D25R75=length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="D25R75"))]) else D25R75=0
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="D75R25"))])>0) D75R25=length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="D75R25"))]) else D75R25=0
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="D25R25"))])>0) D25R25=length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="D25R25"))]) else D25R25=0
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="AVG"))])>0) AVG=length(spe_sub[spe_sub %in% rownames(subset(data_DR_null,data_DR_null$DR_class=="AVG"))]) else AVG=0
    
    #combine all 
    res <- cbind.data.frame(id,AVG,D75R75,D25R25,D75R25,D25R75)
    
    names(res) <- c('cell','AVG','D75R75','D25R25','D75R25','D25R75')
    return(res)
    
    
  },mc.cores = proc))
  
  
  return(Null_funk_all)
}


#Create random matrice where DR_class where randomly distributed among species (keep our number stable) 

#Mammals
#load(file=file.path(results_dir,"mammals","50km","occ_mammals_list.RData"))
#data_DR_randomize_mammals <- lapply(1:1000, function(i) data.frame(sample(data_DR_mammals$DR_class),row.names = rownames(data_DR_mammals)))
#data_DR_randomize_mammals <- lapply(data_DR_randomize_mammals, setNames, "DR_class")
#SES_funk_mammals <- lapply(1:1000,function(i) {Nb.DR_class(ids= names(occ_mammals_list)[1:61618],data_DR_null=data_DR_randomize_mammals[[i]],occ_mammals_list,proc=30)})  
#save(SES_funk_mammals, file = file.path(results_dir,"mammals","50km","SES_funk_mammals.RData"))

#birds
#load(file=file.path(results_dir,"birds","50km","occ_birds_list.RData"))
#data_DR_randomize_birds <- lapply(1:1000, function(i) data.frame(sample(data_DR_birds$DR_class),row.names = rownames(data_DR_birds)))
#data_DR_randomize_birds <- lapply(data_DR_randomize_birds, setNames, "DR_class")
#SES_funk_birds <- lapply(1:1000,function(i) {Nb.DR_class(ids= names(occ_birds_list)[1:61618],data_DR_null=data_DR_randomize_birds[[i]],occ_birds_list,proc=30)})  
#save(SES_funk_birds, file = file.path(results_dir,"birds","50km","SES_funk_birds.RData"))


# Change the first column to rownames and drop this column 
load(file = file.path(results_dir,"mammals","50km","SES_funk_mammals.RData"))
SES_funk <- lapply(SES_funk_mammals, function(df) {df_out <- df[,-1]
rownames(df_out) <- df[[1]]
df_out
})
Null_mean <-  as.data.frame(aaply(laply(SES_funk, as.matrix), c(2, 3), mean))
Null_sd <- as.data.frame(aaply(laply(SES_funk, as.matrix), c(2, 3), sd))

# Compute 
SES_total_mammals <- data.frame(cell=funk_mammals$cell, 
                                D75R75 = (funk_mammals$D75R75 - Null_mean$D75R75)/Null_sd$D75R75,
                                AVG = (funk_mammals$AVG - Null_mean$AVG)/Null_sd$AVG,
                                D25R25 = (funk_mammals$D25R25 - Null_mean$D25R25)/Null_sd$D25R25,
                                D75R25 = (funk_mammals$D75R25 - Null_mean$D75R25)/Null_sd$D75R25,
                                D25R75 = (funk_mammals$D25R75 - Null_mean$D25R75)/Null_sd$D25R75)

save(SES_total_mammals,file = file.path(results_dir,"mammals","50km","SES_total_mammals.RData"))

# Change the first column to rownames and drop this column 
load(file = file.path(results_dir,"birds","50km","SES_funk_birds.RData"))
SES_funk <- lapply(SES_funk_birds, function(df) {df_out <- df[,-1]
rownames(df_out) <- df[[1]]
df_out
})

SES_funk <- data.frame(matrix(unlist(SES_funk), nrow=length(l), byrow=T))


Null_mean <-  as.data.frame(aaply(laply(SES_funk, as.matrix), c(2, 3), mean))
Null_sd <- as.data.frame(aaply(laply(SES_funk, as.matrix), c(2, 3), sd))

# Compute 
SES_total_birds <- data.frame(cell=funk_birds$cell, 
                                D75R75 = (funk_birds$D75R75 - Null_mean$D75R75)/Null_sd$D75R75,
                                AVG = (funk_birds$AVG - Null_mean$AVG)/Null_sd$AVG,
                                D25R25 = (funk_birds$D25R25 - Null_mean$D25R25)/Null_sd$D25R25,
                                D75R25 = (funk_birds$D75R25 - Null_mean$D75R25)/Null_sd$D75R25,
                                D25R75 = (funk_birds$D25R75 - Null_mean$D25R75)/Null_sd$D25R75)

save(SES_total_birds,file = file.path(results_dir,"birds","50km","SES_total_birds.RData"))
