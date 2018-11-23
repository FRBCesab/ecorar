#FUNKRARE

rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NM")

library(funrar)
library(moments)
library(ggplot2)
library(parallel)
library(ade4)
library(dplyr)
library(gridExtra)
library(cluster)
library(rgdal)

#LOAD TRAITS MAPS AND DISTRIB----

  # Load spatial grid for plotting
    map<-readOGR(file.path(data_dir,"ReferenceGrid10Km","gridLand10km.B.shp"))
    #names of each cell
      ID_cell<-rownames(map@data)
    
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
    
  #Check if each species have at least on occurence # TODO a reprendre car il ne check pas, il impose 
    occ_mammals <- occ_mammals[,colSums(occ_mammals)>0]

#----

#COMPUTE FR ----

  ## Mammals
    
    mammalstrait<-mammalstrait[,-1]
    
      ###Format the traits 

        diet <- prep.fuzzy(mammalstrait[,1:10], col.blocks = ncol(mammalstrait[,1:10]), label = "diet")
        
        ForStrat <- mammalstrait$ForStrat %>% as.data.frame()
        row.names(ForStrat) <- row.names(mammalstrait)
        
        Activity <- apply(mammalstrait[,12:14],2,as.numeric) %>% as.data.frame()
        row.names(Activity) <- row.names(mammalstrait)
        Activity=prep.binary(Activity, col.blocks = ncol(Activity), label = "Activity")
        
        bodymass <- log(mammalstrait$BodyMass.Value) %>% as.data.frame()
        row.names(bodymass) <- row.names(mammalstrait)
        bodymass <- bodymass/(max(bodymass, na.rm=T)-min(bodymass, na.rm=T))

    ###Compute the dist matrix
      disTraits_mammals <- dist.ktab(ktab.list.df(list(diet, ForStrat, Activity, bodymass)), c("F","N","B","Q"), scan = FALSE) %>% as.matrix()
      save(disTraits_mammals, file=file.path(results_dir,"mammals/disTraits_mammals.RData"))

    ###Compute funrare indices (note occ_mat are sparse matrices)
      

      FR_mammals <-  funrar(occ_mammals, disTraits_mammals, rel_abund = FALSE) #TODO : pour l'instant cela coince car matrice d'occ trop grande ... 
      
      
      FR_mammals$Ui$species <- as.character(FR_mammals$Ui$species)
      FR_mammals$Ri$species <- as.character(FR_mammals$Ri$species)
      
      Glob_distinc <- distinctiveness_glob(com_dist=disTraits_mammals,abund=NULL)
      colnames(Glob_distinc) <- c("species","Di")
      Glob_distinc$species <- as.character(Glob_distinc$species)
      
      FR_mammals <- list(Ui = FR_mammals$Ui, Di = FR_mammals$Di, Ri = FR_mammals$Ri, GDi=Glob_distinc)
      
      save(FR_mammals, file=file.path(results_dir,"mammals/fr_mammals.RData"))
      

  ##Birds 

    ###Format the traits 

      diet <- prep.fuzzy(birdstrait[,1:10], col.blocks = ncol(birdstrait[,1:10]), label = "diet")
      
      ForStrat <- prep.fuzzy(birdstrait[,11:17], col.blocks = ncol(birdstrait[,11:17]), label = "ForStrat")
      
      bodymass <- log(birdstrait$BodyMass.Value) %>% as.data.frame()
      row.names(bodymass) <- row.names(birdstrait)
      bodymass <- bodymass/(max(bodymass, na.rm=T)-min(bodymass, na.rm=T))
      
      PelagicSpecialist <- as.data.frame(birdstrait$PelagicSpecialist %>% as.character() %>% as.numeric())
      row.names(PelagicSpecialist) <- row.names(birdstrait)
      
      Nocturnal <- as.data.frame(birdstrait$Nocturnal %>% as.character() %>% as.numeric())
      row.names(Nocturnal) <- row.names(birdstrait)

    ###Compute the dist matrix TODO we use scan=FALSE but we probably need to set it to TRUE and chose the best distances ... 

      disTraits_birds <- dist.ktab(ktab.list.df(list(diet, ForStrat,bodymass, PelagicSpecialist, Nocturnal)), c("F","F","Q", "D","D"), scan = FALSE) %>% as.matrix()
      save(disTraits_birds, file=file.path(results_dir,"birds/disTraits_birds.RData"))

    ###Compute funrare indices 

      FR_birds <-  funrar(occ_birds, disTraits_birds, rel_abund = FALSE)
      
      FR_birds$Ui$species <- as.character(FR_birds$Ui$species)
      FR_birds$Ri$species <- as.character(FR_birds$Ri$species)
      
      Glob_distinc <- distinctiveness_glob(com_dist=disTraits_birds,abund=NULL)
      colnames(Glob_distinc) <- c("species","Di")
      Glob_distinc$species <- as.character(Glob_distinc$species)
      
      FR_birds <- list(Ui = FR_birds$Ui, Di = FR_birds$Di, Ri = FR_birds$Ri, GDi=Glob_distinc)
      
      save(FR_birds, file=file.path(results_dir,"birds/fr_birds.RData"))

      
  ## Combine all indices 

    load(file=file.path(results_dir,"mammals/fr_mammals.RData"))
    load(file=file.path(results_dir,"birds/fr_birds.RData"))
    
    comp.fr.all <- function(datatraits,FR_raw,occmat)    
    {
      # datatraits <- mammalstrait
      # FR_raw <- FR_mammals
      # occmat <-occ_mammals 
      
      #calcul the total occurences 
      
      totocc <- data.frame(apply(occ_mammals,2,sum))
      colnames(totocc) <- "occtot"
      
      
      #Compute the FRs 
      
      FR_data <- merge(FR_raw$Ui,FR_raw$Ri)
      FR_data <- merge(FR_data,FR_raw$GDi)
      
      FR_data <- mutate(FR_data, Uin = (Ui-min(Ui)) / max(Ui-min(Ui)),Din = (Di-min(Di)) / max(Di-min(Di)),Rin = (Ri-min(Ri)) / max(Ri-min(Ri)))
      
      FR_data <- mutate(FR_data,FRU=(Uin+Rin)/2,FRD_A=(Din+Rin)/2,FRD_G=(Din*Rin)/2)
      
      #TDsp and rownames
      FR_data$species <- as.character(FR_data$species)
      rownames(FR_data) <- FR_data$species
      
      #Merge with toocc
      
      FR_data <- merge(FR_data,totocc,by=0)
      rownames(FR_data) <- FR_data$species
      FR_data <- FR_data[,-1]
      
      
      # 90% quantile
      Q90_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.1))[10])
      Q90_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.1))[10])
      Q10_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.1))[2])
      Q90_FRD_A <- as.numeric(quantile(FR_data$FRD_A,probs = seq(0, 1, 0.1))[10])
      Q90_FRD_G <- as.numeric(quantile(FR_data$FRD_G,probs = seq(0, 1, 0.1))[10])
      
      # 75% quantile
      Q75_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.25))[4])  
      Q75_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.25))[4])  
      
      # 25% quantile
      Q25_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.25))[2])  
      Q25_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.25))[2])  
      
      Q <- data.frame(Q90_D=Q90_D,Q10_D=Q10_D,Q90_R=Q90_R,Q90_FRD_A=Q90_FRD_A,Q90_FRD_G=Q90_FRD_G,
                      Q75_D=Q75_D,Q75_R=Q75_R,Q25_D=Q25_D,Q25_R=Q25_R)
      
      list(Di=FR_raw$Di,FR=FR_data,Q=Q)
      
    }
    
    FR_birds_all <- comp.fr.all(datatraits=birdstrait,FR_raw=FR_birds)
    save(FR_birds_all, file=file.path(results_dir,"birds/FR_birds_all.RData"))
    
    FR_mammals_all <- comp.fr.all(datatraits=mammalstrait,FR_raw=FR_mammals)
    save(FR_mammals_all, file=file.path(results_dir,"mammals/FR_mammals_all.RData"))

#----
    

#LOAD FR_ALL ---- 
load(file=file.path(results_dir,"birds/FR_birds_all.RData"))
load(file=file.path(results_dir,"mammals/FR_mammals_all.RData"))

#----

#COMPUTE F_SEB ----
source(file.path(script_dir,"multidimFD.R"))
source(file.path(script_dir,"quality_funct_space.R"))
nbdim=4

##Birds
###Compute the qual_funct_space
qual_funct_space_birds<-quality_funct_space(birdstrait, traits_weights=NULL, nbdim=nbdim, metric="Gower", dendro=FALSE,plot=NA)
save(qual_funct_space_birds, file=file.path(results_dir,"birds/qual_funct_space_birds.RData"))
coord_trait_birds<-qual_funct_space_birds$details_funct_space$mat_coord[,1:nbdim]
ggplot(as.data.frame(coord_trait_birds), aes(PC1,PC2)) +geom_point() +ggtitle("Birds")
ggplot(as.data.frame(coord_trait_birds), aes(PC2,PC3)) +geom_point()

###Compute the F_seb metrics 
load(file=file.path(results_dir,"birds/qual_funct_space_birds.RData"))
coord_trait_birds<-qual_funct_space_birds$details_funct_space$mat_coord[,1:nbdim]
F_Seb_birds <- multidimFD(coord_trait_birds, occ_birds, check_species_pool=TRUE, verb=TRUE,nm_asb_plot=NULL)
save(F_Seb_birds, file=file.path(results_dir,"birds/F_Seb_birds.RData"))

#mammals
###Compute the qual_funct_space
qual_funct_space_mammals<-quality_funct_space(mammalstrait, traits_weights=NULL, nbdim=nbdim, metric="Gower", dendro=FALSE,plot=NA)
save(qual_funct_space_mammals, file=file.path(results_dir,"mammals/qual_funct_space_mammals.RData"))
coord_trait_mammals<-qual_funct_space_mammals$details_funct_space$mat_coord[,1:nbdim]
ggplot(as.data.frame(coord_trait_mammals), aes(PC1,PC2)) +geom_point()+ggtitle("mammals")

###Compute the F_seb metrics 
load(file=file.path(results_dir,"mammals/qual_funct_space_mammals.RData"))
coord_trait_mammals<-qual_funct_space_mammals$details_funct_space$mat_coord[,1:nbdim]
F_Seb_mammals <- multidimFD(coord_trait_mammals, occ_mammals, check_species_pool=TRUE, verb=TRUE,nm_asb_plot=NULL)
save(F_Seb_mammals, file=file.path(results_dir,"mammals/F_Seb_mammals.RData"))
#----

#LOAD F_SEB----

load(file.path(results_dir,"mammals/F_Seb_mammals.RData"))
load(file.path(results_dir,"birds/F_Seb_birds.RData"))
#----

#COMPUTE FINAL DATAFRAME---- 

##Generate the subset data 

sub.data <- function(ids,proc,occ_mat,FR_data){
  
  # proc <- 50
  # occ_mat <- occ_birds
  # ids <- rownames(occ_mat)
  # FR_data=FR_birds_all
  
  subD90 <- mclapply(ids,function(id) {    
    spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
    spe_sub[FR_data$FR[spe_sub,"Din"]>FR_data$Q$Q90_D]
    
  },mc.cores = proc)
  names(subD90) <- ids
  
  subR90 <- mclapply(ids,function(id) {    
    spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
    spe_sub[FR_data$FR[spe_sub,"Rin"]>FR_data$Q$Q90_R]
  },mc.cores = proc)
  names(subR90) <- ids
  
  subFRD90_A <- mclapply(ids,function(id) {
    spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
    spe_sub[FR_data$FR[spe_sub,"FRD_A"]>FR_data$Q$Q90_FRD_A]
  },mc.cores = proc)
  names(subFRD90_A) <- ids
  
  subFRD90_G <- mclapply(ids,function(id) {
    spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
    spe_sub[FR_data$FR[spe_sub,"FRD_G"]>FR_data$Q$Q90_FRD_G]
  },mc.cores = proc)
  names(subFRD90_G) <- ids
  
  subDR75 <- mclapply(ids,function(id) {
    spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
    spe_sub_D <- spe_sub[FR_data$FR[spe_sub,"Din"]>FR_data$Q$Q75_D]
    spe_sub_D[FR_data$FR[spe_sub_D,"Rin"]>FR_data$Q$Q75_R]
    
  },mc.cores = proc)
  names(subDR75) <- ids
  
  subDR25 <- mclapply(ids,function(id) {  
    spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
    spe_sub_D <- spe_sub[FR_data$FR[spe_sub,"Din"]<FR_data$Q$Q25_D]
    spe_sub_D[FR_data$FR[spe_sub_D,"Rin"]<FR_data$Q$Q25_R]
  },mc.cores = proc)
  names(subDR25) <- ids
  
  subD75R1 <- mclapply(ids,function(id) {    
    spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
    spe_sub_D <- spe_sub[FR_data$FR[spe_sub,"Din"]>FR_data$Q$Q75_D]
    spe_sub_D[FR_data$FR[spe_sub_D,"Rin"]>0.99]
  },mc.cores = proc)
  names(subD75R1) <- ids
  
  all <- list(subD90,subR90,subFRD90_A,subFRD90_G,subDR75,subDR25,subD75R1)
  names(all) <- c("subD90","subR90","subFRD90_A","subFRD90_G","subDR75","subDR25","subD75R1")
  return(all)
}

sub_birds <- sub.data(ids=rownames(occ_birds),proc=50,occ_mat=occ_birds,FR_data=FR_birds_all)
save(sub_birds, file=file.path(results_dir,"birds/sub_birds.RData"))

sub_mammals <- sub.data(ids=rownames(occ_mammals),proc=50,occ_mat=occ_mammals,FR_data=FR_mammals_all)
save(sub_mammals, file=file.path(results_dir,"mammals/sub_mammals.RData"))

load(file=file.path(results_dir,"birds/sub_birds.RData"))
load(file=file.path(results_dir,"mammals/sub_mammals.RData"))

##Generate main results

final.results <- function(ids,proc,occ_mat,FR_data,F_Seb_data,mat_neigh,seedrand,sub_data){
  
  # proc <- 50
  # occ_mat <- occ_birds
  # ids <- rownames(occ_mat)
  # FR_data <- FR_birds_all
  # F_Seb_data <- F_Seb_birds
  # sub_data <- sub_birds
  
  funk_all <-do.call(rbind,mclapply(ids,function(id) {    #cat("id:",i,"\n")
    
    #id <- ids[1]
    
    TD_sp=sum(occ_mat[id,])
    
    FEve=as.numeric(F_Seb_data[rownames(F_Seb_data)==id,]["FEve"])
    FDis=as.numeric(F_Seb_data[rownames(F_Seb_data)==id,]["FDis"])
    FOri=as.numeric(F_Seb_data[rownames(F_Seb_data)==id,]["FOri"])
    
    FDiv=as.numeric(F_Seb_data[rownames(F_Seb_data)==id,]["FDiv"])
    FRic=as.numeric(F_Seb_data[rownames(F_Seb_data)==id,]["FRic"])
    FSpe=as.numeric(F_Seb_data[rownames(F_Seb_data)==id,]["FSpe"])
    
    #Subset the species present in the site i 
    spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
    
    #Compute the FRD_mean
    
    FRD_A_mean=mean(FR_data$FR[spe_sub,"FRD_A"],na.rm=TRUE)
    FRD_A_skw=skewness(FR_data$FR[spe_sub,"FRD_A"],na.rm=TRUE)
    FRD_G_mean=mean(FR_data$FR[spe_sub,"FRD_G"],na.rm=TRUE)
    FRD_G_skw=skewness(FR_data$FR[spe_sub,"FRD_G"],na.rm=TRUE)
    
    #Number of species within a sites with FR values above the 90% quantile of the whole distribution 
    
    if (length(sub_data$subD90[[id]])>0) D90=length(sub_data$subD90[[id]]) else D90=0
    if (length(sub_data$subR90[[id]])>0) R90=length(sub_data$subR90[[id]]) else R90=0
    if (length(sub_data$subFRD90_A[[id]])>0) FRD90_A=length(sub_data$subFRD90_A[[id]]) else FRD90_A=0
    if (length(sub_data$subFRD90_G[[id]])>0) FRD90_G=length(sub_data$subFRD90_G[[id]]) else FRD90_G=0
    if (length(sub_data$subDR75[[id]])>0) DR75=length(sub_data$subDR75[[id]]) else DR75=0
    if (length(sub_data$subDR25[[id]])>0) DR25=length(sub_data$subDR25[[id]]) else DR25=0
    if (length(sub_data$subD75R1[[id]])>0) D75R1=length(sub_data$subD75R1[[id]]) else D75R1=0
    
    #combine all 
    
    res <- cbind.data.frame(id,TD_sp,FEve,FDis,FOri,FDiv,FRic,FSpe,FRD_A_mean, FRD_A_skw,FRD_G_mean, FRD_G_skw, 
                            D90, R90,FRD90_A, FRD90_G, DR75,DR25,D75R1)
    
    names(res) <- c('poly','TD_sp','FEve','FDis','FOri','FDiv','FRic','FSpe', 'FRD_A_mean', 'FRD_A_skw','FRD_G_mean', 'FRD_G_skw', 
                    'D90','R90','FRD90_A','FRD90_G', 'DR75','DR25','D75R1')
    return(res)
    
  },mc.cores = proc))
  
  return(funk_all)
}

funk_birds <- final.results(ids=rownames(occ_birds),proc=50,occ_mat=occ_birds,FR_data=FR_birds_all,
                            F_Seb_data=F_Seb_birds,mat_neigh=Mat_neighbour_birds,seedrand=1871,sub_data=sub_birds)
save(funk_birds, file=file.path(results_dir,"birds/funk_birds.RData"))

funk_mammals <- final.results(ids=rownames(occ_mammals),proc=50,occ_mat=occ_mammals,FR_data=FR_mammals_all,
                             F_Seb_data=F_Seb_mammals,mat_neigh=Mat_neighbour_mammals,seedrand=1871,sub_data=sub_mammals)
save(funk_mammals, file=file.path(results_dir,"mammals/funk_mammals.RData"))

#----


#VISU ----

load(file=file.path(results_dir,"birds/FR_birds_all.RData"))
load(file=file.path(results_dir,"mammals/FR_mammals_all.RData"))

load(file=file.path(results_dir,"mammals/funk_mammals.RData"))
load(file=file.path(results_dir,"birds/funk_birds.RData"))



##Histograms of species numbers 

ggplot(funk_mammals, aes(TD_sp))+geom_histogram()+ggtitle('mammals')
ggplot(funk_birds, aes(TD_sp))+geom_histogram()+ggtitle('Birds')

##Distributions 

###Birds
a <- ggplot(data=FR_birds_all$FR, aes(x=Din)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_birds_all$Q$Q90_D), colour="red")
b <- ggplot(data=FR_birds_all$FR, aes(x=Rin)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_birds_all$Q$Q90_R), colour="red")
c <- ggplot(data=FR_birds_all$FR, aes(x=FRD_A)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_birds_all$Q$Q90_FRD_A), colour="red")
d <- ggplot(data=FR_birds_all$FR, aes(x=FRD_G)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_birds_all$Q$Q90_FRD_G), colour="red")
grid.arrange(a,b,c,d,ncol=2)

###mammals
a <- ggplot(data=FR_mammals_all$FR, aes(x=Din)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_mammals_all$Q$Q90_D), colour="red")
b <- ggplot(data=FR_mammals_all$FR, aes(x=Rin)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_mammals_all$Q$Q90_R), colour="red")
c <- ggplot(data=FR_mammals_all$FR, aes(x=FRD_A)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_mammals_all$Q$Q90_FRD_A), colour="red")
d <- ggplot(data=FR_mammals_all$FR, aes(x=FRD_G)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_mammals_all$Q$Q90_FRD_G), colour="red")
grid.arrange(a,b,c,d,ncol=2)


##Di vs vs Ri BIRDS 

sub_FR <- FR_birds_all$FR[,c('Rin','Din')]
sub_FR$inv_Rin <- 1-sub_FR$Rin

Q90R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.1))[2])
Q90D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.1))[10])
Q75R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.25))[2])
Q25R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.25))[4])
Q75D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.25))[4])
Q25D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.25))[2])

ggplot(sub_FR, aes(x=inv_Rin,y=Din))+geom_point() + 
  labs(x = "1-Rin",y="Din") +ggtitle('Birds')

ggplot(sub_FR, aes(x=inv_Rin,y=Din))+geom_point() + 
  scale_x_continuous(trans='log10') + 
  geom_vline(aes(xintercept = Q90R), colour="red") + 
  geom_hline(aes(yintercept = Q90D), colour="red") + 
  labs(x = "1-Rin",y="Din")

ggplot(sub_FR, aes(x=inv_Rin,y=Din))+geom_point() + 
  scale_x_continuous(trans='log10') + 
  geom_vline(aes(xintercept = Q75R), colour="red") + 
  geom_hline(aes(yintercept = Q75D), colour="red") + 
  geom_vline(aes(xintercept = Q25R), colour="red") + 
  geom_hline(aes(yintercept = Q25D), colour="red") + 
  labs(x = "1-Rin",y="Din")


##Di vs vs Ri MAMMALS 

sub_FR <- FR_mammals_all$FR[,c('Rin','Din')]
sub_FR$inv_Rin <- 1-sub_FR$Rin

Q90R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.1))[2])
Q90D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.1))[10])
Q75R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.25))[2])
Q25R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.25))[4])
Q75D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.25))[4])
Q25D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.25))[2])



ggplot(sub_FR, aes(x=inv_Rin,y=Din))+geom_point() + 
  labs(x = "1-Rin",y="Din") +ggtitle('Mammals')

ggplot(sub_FR, aes(x=inv_Rin,y=Din))+geom_point() + 
  scale_x_continuous(trans='log10') + 
  geom_vline(aes(xintercept = Q90R), colour="red") + 
  geom_hline(aes(yintercept = Q90D), colour="red") + 
  labs(x = "1-Rin",y="Din")

ggplot(sub_FR, aes(x=inv_Rin,y=Din))+geom_point() + 
  scale_x_continuous(trans='log10') + 
  geom_vline(aes(xintercept = Q75R), colour="red") + 
  geom_hline(aes(yintercept = Q75D), colour="red") + 
  geom_vline(aes(xintercept = Q25R), colour="red") + 
  geom_hline(aes(yintercept = Q25D), colour="red") + 
  labs(x = "1-Rin",y="Din")

##Di Local vs regional 

###Birds Mean 
LDi_temp <- apply(FR_birds_all$Di,2,mean,na.rm=TRUE)
LDi <- data.frame(species=names(LDi_temp),LDi=as.numeric(LDi_temp))

GDi <- data.frame(species=FR_birds_all$FR$species,GDi=FR_birds_all$FR$Di)

Di_all_birds <- merge(LDi,GDi)

a <- ggplot(Di_all_birds, aes(x=GDi,y=LDi))+geom_point() +ggtitle('Birds')

###Birds Var  
LDi_temp <- apply(FR_birds_all$Di,2,sd,na.rm=TRUE)
LDi_sd <- data.frame(species=names(LDi_temp),LDi_sd=as.numeric(LDi_temp))

GDi <- data.frame(species=FR_birds_all$FR$species,GDi=FR_birds_all$FR$Di)

Di_all_birds <- merge(LDi_sd,GDi)

b <- ggplot(Di_all_birds, aes(x=GDi,y=LDi_sd))+geom_point() +ggtitle('Birds')


###mammals Mean 
LDi_temp <- apply(FR_mammals_all$Di,2,mean,na.rm=TRUE)
LDi <- data.frame(species=names(LDi_temp),LDi=as.numeric(LDi_temp))

GDi <- data.frame(species=FR_mammals_all$FR$species,GDi=FR_mammals_all$FR$Di)

Di_all_mammals <- merge(LDi,GDi)

c <- ggplot(Di_all_mammals, aes(x=GDi,y=LDi))+geom_point() +ggtitle('Mammals')

###mammals Var  
LDi_temp <- apply(FR_mammals_all$Di,2,sd,na.rm=TRUE)
LDi_sd <- data.frame(species=names(LDi_temp),LDi_sd=as.numeric(LDi_temp))

GDi <- data.frame(species=FR_mammals_all$FR$species,GDi=FR_mammals_all$FR$Di)

Di_all_mammals <- merge(LDi_sd,GDi)

d <- ggplot(Di_all_mammals, aes(x=GDi,y=LDi_sd))+geom_point() +ggtitle('Mammals')


grid.arrange(a,b,c,d,ncol=2)


##Hist all 

funk_birds$D90
funk_birds$FRD90_A
funk_birds$DR75
funk_birds$DR25

d90 <- data.frame(D90=funk_birds[funk_birds$D90>0,'D90'])
frd90 <- data.frame(FRD90_A=funk_birds[funk_birds$FRD90_A>0,'FRD90_A'])
dr7590 <- data.frame(DR75=funk_birds[funk_birds$DR75>0,'DR75'])
dr2590 <- data.frame(DR25=funk_birds[funk_birds$DR25>0,'DR25'])

a <- ggplot(d90, aes(D90))+geom_histogram(bins=20)
b <- ggplot(frd90, aes(FRD90_A))+geom_histogram(bins=20)
c <- ggplot(dr7590, aes(DR75))+geom_histogram(bins=20)
d <- ggplot(dr2590, aes(DR25))+geom_histogram(bins=20)
grid.arrange(a,b,c,d,ncol=2)

##Corelations

varsub <- c('TD_sp','D90','FRD90_A','FRD90_G','DR75','DR25')

ellip(data=funk_birds,varsub=varsub,plotpdf=FALSE,ordvar='TD_sp',taxa="Birds")
ellip(data=funk_mammals,varsub=varsub,plotpdf=FALSE,ordvar='TD_sp',taxa="mammals")

a <- g.plot(dataplot=funk_birds,x="TD_sp",y="D90",taxa="Birds")
b <- g.plot(dataplot=funk_birds,x="TD_sp",y="FRD90_A",taxa="Birds")
c <- g.plot(dataplot=funk_birds,x="TD_sp",y="DR75",taxa="Birds")
d <- g.plot(dataplot=funk_birds,x="TD_sp",y="DR25",taxa="Birds")
grid.arrange(a,b,c,d,ncol=2)





#Histograms of Ri 

Rinv <- data.frame(Rinv=1-FR_birds_all$FR$Rin)

Q10_Rinv <- as.numeric(quantile(Rinv$Rinv,probs = seq(0, 1, 0.1))[2])

ggplot(Rinv, aes(Rinv))+geom_histogram(bins=20) + scale_x_log10() + geom_vline(aes(xintercept = Q10_Rinv), colour="black")




#Histograms of Di 

a <- ggplot(FR_mammals_all$FR, aes(Din))+geom_histogram(bins=20)+ggtitle('mammals')
b <- ggplot(FR_birds_all$FR, aes(Din))+geom_histogram(bins=20)+ggtitle('Birds')
grid.arrange(a,b,ncol=2)

#Histograms of FRD 

a <- ggplot(FR_mammals_all$FR, aes(Din))+geom_histogram(bins=20)
b <- ggplot(FR_mammals_all$FR, aes(Rin))+geom_histogram(bins=20)
c <- ggplot(FR_mammals_all$FR, aes(FRD_A))+geom_histogram(bins=20)
d <- ggplot(FR_mammals_all$FR, aes(FRD_G))+geom_histogram(bins=20)
grid.arrange(a,b,c,d,ncol=4)

a <- ggplot(FR_birds_all$FR, aes(Din))+geom_histogram(bins=20)
b <- ggplot(FR_birds_all$FR, aes(Rin))+geom_histogram(bins=20)
c <- ggplot(FR_birds_all$FR, aes(FRD_A))+geom_histogram(bins=20)
d <- ggplot(FR_birds_all$FR, aes(FRD_G))+geom_histogram(bins=20)
grid.arrange(a,b,c,d,ncol=4)


#Correlations ... and SES models 

a <- g.plot(dataplot=funk_birds,x="TD_sp",y="D90")
b <- g.plot(dataplot=funk_birds,x="TD_sp",y="D90_SES")
c <- g.plot(dataplot=funk_birds,x="TD_sp",y="FRD90")
d <- g.plot(dataplot=funk_birds,x="TD_sp",y="FRD90_SES")
e <- g.plot(dataplot=funk_birds,x="TD_sp",y="DR80")
f <- g.plot(dataplot=funk_birds,x="TD_sp",y="DR80_SES")
g <- g.plot(dataplot=funk_birds,x="TD_sp",y="R90")
h <- g.plot(dataplot=funk_birds,x="TD_sp",y="R90_SES")

grid.arrange(a,g,c,e,b,h,d,f,ncol=4)


FR_birds_all$FR$
  
  min(funk_birds$UiSkw)
max(funk_birds$UiSkw)

fb_pos <- funk_birds$UiSkw[funk_birds$UiSkw>0]
fb_pos <- fb_pos[complete.cases(fb_pos)]


set.seed(2)

flow1 <- funk_birds[funk_birds$UiSkw<0.5,]
id_lowskew <- flow1[order(flow1$TD_sp,decreasing = TRUE),"poly"][1]
spe_sub_lowskew <- names(occ_birds[id_lowskew,][occ_birds[id_lowskew,]>0])
hist(FR_birds$Ui[FR_birds$Ui$species %in% spe_sub_lowskew,2])

flow2 <- funk_birds[funk_birds$UiSkw>2.75,]
id_highskew <- flow2[order(flow2$TD_sp,decreasing = TRUE),"poly"][33]
spe_sub_highskew <- names(occ_birds[id_highskew,][occ_birds[id_highskew,]>0])
Ui_highskew <- FR_birds$Ui[FR_birds$Ui$species %in% spe_sub_highskew,]
Di_highskew <-  as.data.frame(FR_birds$Di[id_highskew,])
colnames(Di_highskew) <- c('Di')
Ri_highskew <- FR_birds$Ri[FR_birds$Ri$species %in% spe_sub_highskew,]

ggplot(data=as.data.frame(Ui_highskew), aes(x=Ui)) + geom_histogram(binwidth=0.01,colour = "black") + geom_vline(aes(xintercept = Q_Ui_birds), colour="black")
ggplot(data=as.data.frame(Di_highskew), aes(x=Di)) + geom_histogram(binwidth=0.01,colour = "black")
ggplot(data=as.data.frame(Ri_highskew), aes(x=Ri)) + geom_histogram(binwidth=0.01,colour = "black")
#----    








