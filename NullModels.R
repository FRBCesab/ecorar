#Null model
rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NL")

FR_mammals

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


##Generate main results
SES.final.results <- function(ids,proc,occ_mat_list,sub_data){
  
  # proc <- 50
  # occ_mat_list <- occ_mammals_list
  # ids <- names(occ_mat_list)
  # sub_data <- sub_mammals
  
  SES_funk_all <-do.call(rbind,mclapply(ids,function(id) {    #cat("id:",i,"\n")
    
    #id <- ids[10000]
    
    #Subset the species present in the site i 
    spe_sub <- unique(occ_mat_list[[id]])
    
    #Number of species within a sites with FR values above quantiles of species and functional distribution 
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D75R75"))])>0) D75R75=length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D75R75"))]) else D75R75=0
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D25R75"))])>0) D25R75=length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D25R75"))]) else D25R75=0
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D75R25"))])>0) D75R25=length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D75R25"))]) else D75R25=0
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D25R25"))])>0) D25R25=length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D25R25"))]) else D25R25=0
    if (length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="AVG"))])>0) AVG=length(spe_sub[spe_sub %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="AVG"))]) else AVG=0
    
    #combine all 
    res <- cbind.data.frame(id,AVG,D75R75,D25R25,D75R25,D25R75)
    
    names(res) <- c('cell','AVG','D75R75','D25R25','D75R25','D25R75')
    return(res)
    
  },mc.cores = proc))
  
  return(funk_all)
}

data_DR_mammals


SES_funk_birds <-
SES_funk_mammals <- 

(obs - mean(null))/sd(null) )