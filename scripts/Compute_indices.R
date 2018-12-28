#FUNKRARE

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

#LOAD TRAITS MAPS AND DISTRIB----
reso="50km"
  # Load traits and distrib 
    load(file=file.path(results_dir,"mammals",reso,"mammalsID.RData"))
    load(file=file.path(results_dir,"mammals",reso,"mammalstrait.RData"))
    load(file=file.path(results_dir,"mammals",reso,"occ_mammals_list.RData"))
#----
    mammalsID<-mammalsID[mammalsID$ID %in% rownames(mammalstrait),]
    
#COMPUTE FR ----

## Mammals
 
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
      save(disTraits_mammals, file=file.path(results_dir,"mammals",reso,"disTraits_mammals.RData"))

    ###Compute funrare indices (note occ_mat are sparse matrices)
      
      load(file=file.path(results_dir,"mammals",reso,"disTraits_mammals.RData"))
      
     #matrice to big, build hypothetical community where all species are presents. Allow to compute Ui & Di for each species
      Sim_commu <- matrix(1,1,ncol(disTraits_mammals))
      colnames(Sim_commu) <- colnames(disTraits_mammals)
      
      #Compute Ui (Global Uniqueness)
        Ui<-uniqueness(Sim_commu,disTraits_mammals)
        rownames(Ui)<-Ui[,1]
      #Compute Di (Global distinctiveness)
        Di<-t(distinctiveness(Sim_commu,disTraits_mammals))
        colnames(Di)<-"Di"
     
      #matrix is to big to compute Ri + restrictedness function need at least 2 species to be compute.
        Ri<-data.frame(table(unlist(occ_mammals_list))/length(occ_mammals_list))
        rownames(Ri)<-Ri[,1]
        Ri<-Ri[,-1, drop = FALSE]
        Ri<-1-Ri
        colnames(Ri)<-"Ri"

        #Create the FR_data frame 
        FR_data <- merge(Ui,Di, by="row.names")
        rownames(FR_data) <- FR_data[,1]
        FR_data <- FR_data[,-c(1,2)]
        FR_data <- merge(FR_data,Ri, by="row.names")
        rownames(FR_data) <- FR_data[,1]
        FR_data <- FR_data[,-1]
          
        FR_data$Uin<-(FR_data$Ui-min(FR_data$Ui)) / max(FR_data$Ui-min(FR_data$Ui))
        FR_data$Din<-(FR_data$Di-min(FR_data$Di)) / max(FR_data$Di-min(FR_data$Di))
        FR_data$Rin<-(FR_data$Ri-min(FR_data$Ri)) / max(FR_data$Ri-min(FR_data$Ri))
  
        # 90% quantile
        Q90_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.1))[10])
        Q90_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.1))[10])
        Q10_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.1))[2])
  
        # 75% quantile
        Q75_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.25))[4])  
        Q75_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.25))[4])  
        
        # 25% quantile
        Q25_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.25))[2])  
        Q25_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.25))[2])  
        
        Q <- data.frame(Q90_D=Q90_D,Q10_D=Q10_D,Q90_R=Q90_R,
                        Q75_D=Q75_D,Q75_R=Q75_R,Q25_D=Q25_D,Q25_R=Q25_R)
        
        FR_mammals <- list(FR=FR_data,Q=Q)
        save(FR_mammals, file=file.path(results_dir,"mammals",reso,"FR_mammals.RData"))
      

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

      ##Compute funrare indices (note occ_mat are sparse matrices)
      
      load(file=file.path(results_dir,"birds/disTraits_birds.RData"))
      
      #matrice to big, build hypothetical community where all species are presents. Allow to compute Ui & Di for each species
      Sim_commu <- matrix(1,1,dim(occ_birds)[2])
      colnames(Sim_commu) <- colnames(occ_birds)
      
      #Compute Ui (Global Uniqueness)
      Ui<-uniqueness(Sim_commu,disTraits_birds)
      
      #Compute Di (Global distinctiveness)
      Di<-t(distinctiveness(Sim_commu,disTraits_birds))
      
      # matrix is to big to compute Ri + restrictedness function need at least 2 species to be compute.
      Ri<-data.frame(1-(colSums(occ_birds)/dim(occ_birds)[1]))
      colnames(Ri)<-"Ri"
      

      #Create the FR_data frame 
      
      FR_data <- data.frame(Ui,Di,Ri)

      FR_data <- mutate(FR_data, Uin = (Ui-min(Ui)) / max(Ui-min(Ui)),Din = (Di-min(Di)) / max(Di-min(Di)),Rin = (Ri-min(Ri)) / max(Ri-min(Ri)))
      
      
      # 90% quantile
      Q90_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.1))[10])
      Q90_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.1))[10])
      Q10_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.1))[2])
      
      # 75% quantile
      Q75_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.25))[4])  
      Q75_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.25))[4])  
      
      # 25% quantile
      Q25_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.25))[2])  
      Q25_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.25))[2])  
      
      Q <- data.frame(Q90_D=Q90_D,Q10_D=Q10_D,Q90_R=Q90_R,
                      Q75_D=Q75_D,Q75_R=Q75_R,Q25_D=Q25_D,Q25_R=Q25_R)
      
      FR_birds <- list(FR=FR_data,Q=Q)
      save(FR_birds, file=file.path(results_dir,"birds/FR_birds.RData"))
#----


      
#COMPUTE FINAL DATAFRAME---- 
      
    load(file=file.path(results_dir,"mammals/FR_mammals.RData"))
      
       ##Generate the subset data 
      
      sub.data <- function(ids,proc,occ_mat_list,FR_data){
        
        #proc <- 3
        #occ_mat_list <- occ_mammals_list
        #ids <- names(occ_mammals_list)
        #FR_data=FR_mammals
        
        # subD90 <- mclapply(ids,function(id) {    
        #   spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
        #   spe_sub[FR_data$FR[spe_sub,"Din"]>FR_data$Q$Q90_D]
        #   
        # },mc.cores = proc)
        # names(subD90) <- ids
        # 
        # subR90 <- mclapply(ids,function(id) {    
        #   spe_sub <- names(occ_mat[id,][occ_mat[id,]>0])
        #   spe_sub[FR_data$FR[spe_sub,"Rin"]>FR_data$Q$Q90_R]
        # },mc.cores = proc)
        # names(subR90) <- ids

        subD75R75 <- mclapply(ids,function(id) {
          spe_sub <- occ_mat_list[[id]]
          spe_sub_D <- spe_sub[FR_data$FR[spe_sub,"Din"]>FR_data$Q$Q75_D]
          spe_sub_D[FR_data$FR[spe_sub_D,"Rin"]>FR_data$Q$Q75_R]
          
        },mc.cores = proc)
        names(subD75R75) <- ids
        
        subD25R25 <- mclapply(ids,function(id) {  
          spe_sub <-  occ_mat_list[[id]]
          spe_sub_D <- spe_sub[FR_data$FR[spe_sub,"Din"]<FR_data$Q$Q25_D]
          spe_sub_D[FR_data$FR[spe_sub_D,"Rin"]<FR_data$Q$Q25_R]
        },mc.cores = proc)
        names(subD25R25) <- ids
        
        subD75R25 <- mclapply(ids,function(id) {
          spe_sub <-  occ_mat_list[[id]]
          spe_sub_D <- spe_sub[FR_data$FR[spe_sub,"Din"]>FR_data$Q$Q75_D]
          spe_sub_D[FR_data$FR[spe_sub_D,"Rin"]<FR_data$Q$Q25_R]
        },mc.cores = proc)
        names(subD75R25) <- ids
        
        subD25R75 <- mclapply(ids,function(id) {  
          spe_sub <-  occ_mat_list[[id]]
          spe_sub_D <- spe_sub[FR_data$FR[spe_sub,"Din"]<FR_data$Q$Q25_D]
          spe_sub_D[FR_data$FR[spe_sub_D,"Rin"]>FR_data$Q$Q75_R]
        },mc.cores = proc)
        names(subD25R75) <- ids
        
        #TODO : AJOUTER LES AVERGAE
        subAVG <- mclapply(ids,function(id) {  
          spe_sub <-  occ_mat_list[[id]]
          spe_sub_D <- spe_sub[(FR_data$FR[spe_sub,"Din"]>FR_data$Q$Q25_D)&(FR_data$FR[spe_sub,"Din"]<FR_data$Q$Q75_D)]
          spe_sub_D[(FR_data$FR[spe_sub_D,"Rin"]>FR_data$Q$Q25_R) & (FR_data$FR[spe_sub_D,"Rin"]<FR_data$Q$Q75_R)]
        },mc.cores = proc)
        names(subAVG) <- ids
        
        all <- list(subAVG,subD75R75,subD25R25,subD75R25,subD25R75)
        names(all) <- c("subAVG","subD75R75","subD25R25","subD75R25","subD25R75")
        return(all)
      }

sub_mammals <- sub.data(ids=names(occ_mammals_list),proc=3,occ_mat_list=occ_mammals_list,FR_data=FR_mammals)
save(sub_mammals, file=file.path(results_dir,"mammals",reso,"sub_mammals.RData"))


sub_birds <- sub.data(ids=rownames(occ_birds),proc=50,occ_mat=occ_birds,FR_data=FR_birds_all)
save(sub_birds, file=file.path(results_dir,"birds/sub_birds.RData"))


load(file=file.path(results_dir,"birds/sub_birds.RData"))
load(file=file.path(results_dir,"mammals/sub_mammals.RData"))

##Generate main results
final.results <- function(ids,proc,occ_mat_list,sub_data){
  
  # proc <- 50
  # occ_mat_list <- occ_mammals_list
  # ids <- names(occ_mat_list)
  # sub_data <- sub_mammals
  
  funk_all <-do.call(rbind,mclapply(ids,function(id) {    #cat("id:",i,"\n")
    
    #id <- ids[1]
    
    TD_sp=length(unique(occ_mat_list[[id]]))
    
    #Subset the species present in the site i 
    spe_sub <- unique(occ_mat_list[[id]])
    
    #Number of species within a sites with FR values above quantiles of species and functional distribution 
    if (length(sub_data$subD75R75[[id]])>0) D75R75=length(sub_data$subD75R75[[id]]) else D75R75=0
    if (length(sub_data$subD25R25[[id]])>0) D25R25=length(sub_data$subD25R25[[id]]) else D25R25=0
    if (length(sub_data$subD75R25[[id]])>0) D75R25=length(sub_data$subD75R25[[id]]) else D75R25=0
    if (length(sub_data$subD25R75[[id]])>0) D25R75=length(sub_data$subD25R75[[id]]) else D25R75=0
    if (length(sub_data$subAVG[[id]])>0) AVG=length(sub_data$subAVG[[id]]) else AVG=0
   
     #combine all 
    res <- cbind.data.frame(id,TD_sp,AVG,D75R75,D25R25,D75R25,D25R75)
    
    names(res) <- c('cell','TD_sp','AVG','D75R75','D25R25','D75R25','D25R75')
    return(res)
    
  },mc.cores = proc))
  
  return(funk_all)
}

funk_mammals <- final.results(ids=names(occ_mammals_list),proc=3,occ_mat_list=occ_mammals_list,sub_data=sub_mammals)
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








