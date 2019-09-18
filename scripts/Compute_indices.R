#FUNKRARE

rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=TRUE,who="NM")

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
library(sjstats)


#COMPUTE FR ----

  ## MAMMALS

    ### Load traits and distrib 
        load(file=file.path(results_dir,"mammals/mammalsID.RData"))
        load(file=file.path(results_dir,"mammals","mammalstrait.RData"))
        load(file=file.path(results_dir,"mammals","50km","occ_mammals_list.RData"))
        mammalsID<-mammalsID[mammalsID$ID %in% rownames(mammalstrait),]

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
        save(disTraits_mammals, file=file.path(results_dir,"mammals","disTraits_mammals.RData"))

    ###Compute funrare indices (note occ_mat are sparse matrices)
      
        load(file=file.path(results_dir,"mammals","disTraits_mammals.RData"))
      
        ####matrice to big, build hypothetical community where all species are presents. Allow to compute Ui & Di for each species
          Sim_commu <- matrix(1,1,ncol(disTraits_mammals))
          colnames(Sim_commu) <- colnames(disTraits_mammals)
      
        ####Compute Ui (Global Uniqueness)
          Ui<-uniqueness(str(),disTraits_mammals)
          rownames(Ui)<-Ui[,1]
          
        ####Compute Di (Global distinctiveness)
          Di<-t(distinctiveness(Sim_commu,disTraits_mammals))
          colnames(Di)<-"Di"
          
        ####Compute Di_local (Local distinctiveness)
          
          spnames <- unique(unlist(occ_mammals_list))
          
          proc=45
          Di_locall <- mclapply(1:length(occ_mammals_list),function(i){
            id <- occ_mammals_list[[i]]
            if(!is.na(id[2])==TRUE) {  #compute Di only for communities with 2 species or more 
              com <- matrix(data = 0, nrow = 1, ncol = length(spnames))
              colnames(com) <- spnames
              for (j in 1:length(id)) { com[1,id[j]]=1 }
              Di_local <- t(distinctiveness(com,disTraits_mammals))
              colnames(Di_local)<-names(occ_mammals_list[i])
              Di_local <- as.data.frame(Di_local)
              Di_local
            }
 
          },mc.cores = proc)
          
          Di_locall <- plyr::compact(Di_locall) #remove the NULL
          Di_locall <- do.call(cbind,Di_locall)
          
          save(Di_locall, file=file.path(results_dir,"mammals","50km","Di_locall_mammals.RData"))
          
          load(file=file.path(results_dir,"mammals","50km","Di_locall_mammals.RData"))
          mean_loc_Di <- apply(Di_locall,1,function(x) mean(x, na.rm = T)) 
          mean_loc_Di <- data.frame(mean_loc_Di)
          
          sd_loc_Di <- apply(Di_locall,1,function(x) var(x, na.rm = T)) #using sjstats
          sd_loc_Di <- data.frame(sd_loc_Di)

          Di_loc_glob <- merge(Di,mean_loc_Di,by="row.names",all.x=FALSE)
          rownames(Di_loc_glob) <- Di_loc_glob[,1]
          Di_loc_glob <- Di_loc_glob[,-1]
          Di_loc_glob <- merge(Di_loc_glob,sd_loc_Di,by="row.names",all.x=FALSE)
          colnames(Di_loc_glob) <- c("ID","Di_glob","Di_loc_mean","Di_loc_sd")
          
          Di_loc_glob$Di_loc_cv=Di_loc_glob$Di_loc_sd/Di_loc_glob$Di_loc_mean
          
          
          ggplot(Di_loc_glob, aes(x=Di_glob, y=Di_loc_mean)) + 
            geom_point(size=1) + 
            stat_smooth(method = "lm", formula = y ~ x, size = 1,se=TRUE) + 
            ggtitle("Mammals") + 
            geom_abline(intercept = 0, slope = 1, color="red", 
                          linetype="dashed", size=1)+theme_bw()+
            xlab("Global distinctiveness")+ ylab("Local distinctiveness")
          
          ggplot(Di_loc_glob, aes(x=Di_glob, y=Di_loc_var)) + 
            geom_point(size=1) + 
            stat_smooth(method = "lm", formula = y ~ x, size = 1,se=TRUE) + 
            ggtitle("Mammals") + 
            geom_abline(intercept = 0, slope = 1, color="red", 
                        linetype="dashed", size=1)
          
          ggplot(Di_loc_glob, aes(x=Di_glob, y=log(Di_loc_cv))) + 
            geom_point(size=1) + 
            stat_smooth(method = "lm", formula = y ~ x, size = 1,se=TRUE) + 
            ggtitle("Mammals") + 
            geom_abline(intercept = 0, slope = 1, color="red", 
                        linetype="dashed", size=1)
          
          
          a <- ggplot(Di_loc_glob_mean, aes(x=Di_glob, y=Di_loc)) + 
            stat_smooth(method = "lm", formula = y ~ x, size = 1,se=TRUE) + ylim(0, 1) +
            ggtitle("Mamals") + 
            geom_abline(intercept = 0, slope = 1, color="red", 
                        linetype="dashed", size=1)
          
          
          
          Di_loc_glob_all <- merge(Di,Di_locall,by="row.names",all.x=FALSE)
          
          proc=45
          Di_loc_glob_mamals <- do.call(rbind,mclapply(1:dim(Di_loc_glob_all)[1],function(i) {
            one_sp <- Di_loc_glob_all[i,]
            one_sp <- one_sp[,!is.na(one_sp)]
            dim <- length(one_sp[-c(1,2)])
            loc_glog <- data.frame(rep(as.numeric(one_sp[2]),dim),as.numeric(one_sp[-c(1,2)]))
            colnames(loc_glog) <- c("Di_glob","Di_loc")
            loc_glog
          },mc.cores = proc))
          
          save(Di_loc_glob_mamals, file=file.path(results_dir,"mammals","50km","Di_loc_glob_mamals.RData"))

          plot(Di_loc_glob_mamals$Di_glob,Di_loc_glob_mamals$Di_loc)
          
          
        #####matrix is to big to compute Ri + restrictedness function need at least 2 species to be compute.
          Ri<-data.frame(table(unlist(occ_mammals_list))/length(occ_mammals_list))
          rownames(Ri)<-Ri[,1]
          Ri<-Ri[,-1, drop = FALSE]
          Ri<-1-Ri
          colnames(Ri)<-"Ri"

        #####Create the FR_data frame 
          FR_data <- merge(Ui,Di, by="row.names")
          rownames(FR_data) <- FR_data[,1]
          FR_data <- FR_data[,-c(1,2)]
          FR_data <- merge(FR_data,Ri, by="row.names")
          rownames(FR_data) <- FR_data[,1]
          FR_data <- FR_data[,-1]
          
          FR_data$Uin<-(FR_data$Ui-min(FR_data$Ui)) / max(FR_data$Ui-min(FR_data$Ui))
          FR_data$Din<-(FR_data$Di-min(FR_data$Di)) / max(FR_data$Di-min(FR_data$Di))
          FR_data$Rin<-(FR_data$Ri-min(FR_data$Ri)) / max(FR_data$Ri-min(FR_data$Ri))
  
        ##### 90% quantile
          Q90_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.1))[10])
          Q90_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.1))[10])
          Q10_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.1))[2])
  
        ##### 75% quantile
          Q75_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.25))[4])  
          Q75_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.25))[4])  
        
        ##### 25% quantile
          Q25_D <- as.numeric(quantile(FR_data$Din,probs = seq(0, 1, 0.25))[2])  
          Q25_R <- as.numeric(quantile(FR_data$Rin,probs = seq(0, 1, 0.25))[2])  
        
        Q <- data.frame(Q90_D=Q90_D,Q10_D=Q10_D,Q90_R=Q90_R,
                        Q75_D=Q75_D,Q75_R=Q75_R,Q25_D=Q25_D,Q25_R=Q25_R)
        
        FR_mammals <- list(FR=FR_data,Q=Q)
        save(FR_mammals, file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))
      

  ##BIRDS
        
        # Load traits and distrib 
        load(file=file.path(data_dir,"birds","birdsID.RData"))
        load(file=file.path(results_dir,"birds","birdstrait.RData"))
        load(file=file.path(results_dir,"birds","50km","occ_birds_list.RData"))
        
        #----
     
        #COMPUTE FR ----
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
      save(disTraits_birds, file=file.path(results_dir,"birds","50km","disTraits_birds.RData"))

      load(file=file.path(results_dir,"birds","disTraits_birds.RData"))
      
      ##Compute funrare indices (note occ_mat are sparse matrices)
      #matrice to big, build hypothetical community where all species are presents. Allow to compute Ui & Di for each species
      Sim_commu <- matrix(1,1,ncol(disTraits_birds))
      colnames(Sim_commu) <- colnames(disTraits_birds)
      
      #Compute Ui (Global Uniqueness)
      Ui<-uniqueness(Sim_commu,disTraits_birds)
      rownames(Ui)<-Ui[,1]
      #Compute Di (Global distinctiveness)
      Di<-t(distinctiveness(Sim_commu,disTraits_birds))
      colnames(Di)<-"Di"
      
      ####Compute Di_local (Local distinctiveness)
      
      spnames <- unique(unlist(occ_birds_list))
      
      proc=40
      Di_locall <- mclapply(1:length(occ_birds_list),function(i){
        id <- occ_birds_list[[i]]
        if(!is.na(id[2])==TRUE) {  #compute Di only for communities with 2 species or more 
          com <- matrix(data = 0, nrow = 1, ncol = length(spnames))
          colnames(com) <- spnames
          for (j in 1:length(id)) { com[1,id[j]]=1 }
          Di_local <- t(distinctiveness(com,disTraits_birds))
          colnames(Di_local)<-names(occ_birds_list[i])
          Di_local <- as.data.frame(Di_local)
          Di_local
        }
        
      },mc.cores = proc)
      
      Di_locall <- plyr::compact(Di_locall) #remove the NULL
      Di_locall <- do.call(cbind,Di_locall)
      
      save(Di_locall, file=file.path(results_dir,"birds","50km","Di_locall_birds.RData"))
      
      load(file=file.path(results_dir,"birds","50km","Di_locall_birds.RData"))
      mean_loc_Di <- apply(Di_locall,1,function(x) mean(x, na.rm = T)) 
      mean_loc_Di <- data.frame(mean_loc_Di)
      
      sd_loc_Di <- apply(Di_locall,1,function(x) var(x, na.rm = T)) #using sjstats
      sd_loc_Di <- data.frame(sd_loc_Di)
      
      Di_loc_glob <- merge(Di,mean_loc_Di,by="row.names",all.x=FALSE)
      rownames(Di_loc_glob) <- Di_loc_glob[,1]
      Di_loc_glob <- Di_loc_glob[,-1]
      Di_loc_glob <- merge(Di_loc_glob,sd_loc_Di,by="row.names",all.x=FALSE)
      colnames(Di_loc_glob) <- c("ID","Di_glob","Di_loc_mean","Di_loc_sd")
      
      Di_loc_glob$Di_loc_cv=Di_loc_glob$Di_loc_sd/Di_loc_glob$Di_loc_mean
      
      
      
      ggplot(Di_loc_glob, aes(x=Di_glob, y=Di_loc_mean)) + 
        geom_point(size=1) + 
        stat_smooth(method = "lm", formula = y ~ x, size = 1,se=TRUE) + 
        ggtitle("Birds") + 
        geom_abline(intercept = 0, slope = 1, color="red", 
                    linetype="dashed", size=1)+theme_bw()+
        xlab("Global distinctiveness")+ ylab("Local distinctiveness")
    
        #Extract outliers species 
      Di_loc_glob <- Di_loc_glob[Di_loc_glob$ID %in% birdsID$ID,]
      
      
      up <- Di_loc_glob[Di_loc_glob[,"Di_glob"]<0.55,]
      up <- up[up[,"Di_loc"]>0.7,]
      upspecies <- birdsID[birdsID$ID %in% up$ID,]
      
      down <- Di_loc_glob[Di_loc_glob[,"Di_glob"]>0.75,]
      down <- down[down[,"Di_loc"]<0.5,]
      downspecies <- birdsID[birdsID$ID %in% down$ID,]
      
      
      ids <- upspecies$ID
      subUP <- mclapply(ids,function(id) {  
        spe_up <-  names(unlist(lapply(occ_birds_list, function(x) x[x==id])))
      },mc.cores = 4)
      names(subUP) <- ids
      
      
      ids <- downspecies$ID
      subdown <- mclapply(ids,function(id) {  
        spe_down <-  names(unlist(lapply(occ_birds_list, function(x) x[x==id])))
      },mc.cores = 4)
      names(subdown) <- ids
      
      #matrix is to big to compute Ri + restrictedness function need at least 2 species to be compute.
      Ri<-data.frame(table(unlist(occ_birds_list))/length(occ_birds_list))
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
      
      FR_birds <- list(FR=FR_data,Q=Q)
      save(FR_birds, file=file.path(results_dir,"birds","50km","FR_birds.RData"))
      
#----


      
#COMPUTE FINAL DATAFRAME---- 

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

load(file=file.path(results_dir,"mammals","50km","sub_mammals.RData"))
load(file=file.path(results_dir,"birds","50km","sub_birds.RData"))

sub_mammals <- sub.data(ids=names(occ_mammals_list),proc=3,occ_mat_list=occ_mammals_list,FR_data=FR_mammals)
save(sub_mammals, file=file.path(results_dir,"mammals","50km","sub_mammals.RData"))

sub_birds <- sub.data(ids=names(occ_birds_list),proc=3,occ_mat_list=occ_birds_list,FR_data=FR_birds)
save(sub_birds, file=file.path(results_dir,"birds","50km","sub_birds.RData"))

load(file=file.path(results_dir,"birds","50km","sub_birds.RData"))
load(file=file.path(results_dir,"mammals","50km","sub_mammals.RData"))

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
save(funk_mammals, file=file.path(results_dir,"mammals","50km","funk_mammals.RData"))

funk_birds <- final.results(ids=names(occ_birds_list),proc=3,occ_mat_list=occ_birds_list,sub_data=sub_birds)
save(funk_birds, file=file.path(results_dir,"birds","50km","funk_birds.RData"))




#----

load(file=file.path(results_dir,"birds","50km","funk_birds.RData"))
load(file=file.path(results_dir,"mammals","50km","funk_mammals.RData"))

load(file=file.path(results_dir,"birds","50km","FR_birds.RData"))
load(file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))

str(FR_mammals)

#VISU ----
##Histograms of species numbers 

ggplot(funk_mammals, aes(TD_sp))+geom_histogram()+ggtitle('mammals')
ggplot(funk_birds, aes(TD_sp))+geom_histogram()+ggtitle('Birds')

##Distributions 

###Birds
a <- ggplot(data=FR_birds$FR, aes(x=Din)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_birds$Q$Q75_D), colour="red")
b <- ggplot(data=FR_birds$FR, aes(x=Rin)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_birds$Q$Q75_R), colour="red")

grid.arrange(a,b,ncol=2)

###mammals
a <- ggplot(data=FR_mammals$FR, aes(x=Din)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_mammals$Q$Q75_D), colour="red")
b <- ggplot(data=FR_mammals$FR, aes(x=Rin)) + geom_histogram(binwidth=0.05,colour = "black")+ geom_vline(aes(xintercept = FR_mammals$Q$Q75_R), colour="red")
grid.arrange(a,b,ncol=2)


##Di vs vs Ri BIRDS 

sub_FR <- FR_birds$FR[,c('Rin','Din')]
sub_FR$inv_Rin <- 1-sub_FR$Rin

sub_FR <- merge (sub_FR, data_DR_birds, by="row.names")
levels(sub_FR$DR_class) <- c(levels(sub_FR$DR_class), "Average","Rare","Common","Other")
sub_FR$DR_class[sub_FR$DR_class=="AVG"]<-"Average"
sub_FR$DR_class[sub_FR$DR_class=="D75R25"]<-"Other"
sub_FR$DR_class[sub_FR$DR_class=="D75R75"]<-"Rare"
sub_FR$DR_class[sub_FR$DR_class=="D25R75"]<-"Other"
sub_FR$DR_class[sub_FR$DR_class=="D25R25"]<-"Common"
sub_FR$DR_class[is.na(sub_FR$DR_class)]<-"Other"


Q75R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.25))[2])
Q25R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.25))[4])
Q75D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.25))[4])
Q25D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.25))[2])

ggplot(sub_FR, aes(x=inv_Rin,y=Din))+geom_point() + 
  labs(x = "1-Rin",y="Din") +ggtitle('Birds')

ggplot(sub_FR, aes(x=inv_Rin,y=Din,col=DR_class))+geom_point() + 
  scale_x_continuous(trans='log10') + 
  scale_color_manual(values=c("#00AFBB", "#FF4500", "#E7B800", "grey"))+
  geom_vline(aes(xintercept = Q75R), colour="red") + 
  geom_hline(aes(yintercept = Q75D), colour="red") + 
  geom_vline(aes(xintercept = Q25R), colour="red") + 
  geom_hline(aes(yintercept = Q25D), colour="red") + 
  labs(x = "1-Rin",y="Din")+theme_bw()


##Di vs vs Ri MAMMALS 
sub_FR <- FR_mammals$FR[,c('Rin','Din')]
sub_FR$inv_Rin <- 1-sub_FR$Rin

sub_FR <- merge (sub_FR, data_DR_mammals, by="row.names")
levels(sub_FR$DR_class) <- c(levels(sub_FR$DR_class), "Average","Rare","Common","Other")
sub_FR$DR_class[sub_FR$DR_class=="AVG"]<-"Average"
sub_FR$DR_class[sub_FR$DR_class=="D75R25"]<-"Other"
sub_FR$DR_class[sub_FR$DR_class=="D75R75"]<-"Rare"
sub_FR$DR_class[sub_FR$DR_class=="D25R75"]<-"Other"
sub_FR$DR_class[sub_FR$DR_class=="D25R25"]<-"Common"
sub_FR$DR_class[is.na(sub_FR$DR_class)]<-"Other"
Q75R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.25))[2])
Q25R <- as.numeric(quantile(sub_FR$inv_Rin,probs = seq(0, 1, 0.25))[4])
Q75D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.25))[4])
Q25D <- as.numeric(quantile(sub_FR$Din,probs = seq(0, 1, 0.25))[2])

ggplot(sub_FR, aes(x=inv_Rin,y=Din))+geom_point() + 
  labs(x = "1-Rin",y="Din") +ggtitle('Mammals')

ggplot(sub_FR, aes(x=inv_Rin,y=Din,col=DR_class))+geom_point() + 
  scale_x_continuous(trans='log10') + 
  scale_color_manual(values=c("#00AFBB", "#FF4500", "#E7B800", "grey"))+
  geom_vline(aes(xintercept = Q75R), colour="red") + 
  geom_hline(aes(yintercept = Q75D), colour="red") + 
  geom_vline(aes(xintercept = Q25R), colour="red") + 
  geom_hline(aes(yintercept = Q25D), colour="red") + 
  labs(x = "1-Rin",y="Din")+theme_bw()

##Di Local vs regional 
#Compute Di at the local scale


id= mammalsID$ID[2]

  list.2 <- names(unlist(lapply(occ_mammals_list, function(x) subset(x,x==id))))
  list.3 <- occ_mammals_list
  
DiLocal <- function (ids, occ_list, )

  site_sp <- lapply(1:length(occ_mammals_list), function(x) subset(occ_mammals_list[[x]],occ_mammals_list[[x]]==id))
  names(site_sp) <-names(occ_mammals_list)
  site_sp <- c(names(unlist(site_sp)))

 # test <- lapply(site,function(id) {cooccsp<-  occ_mammals_list[[id]]})

  Sim_commu <- matrix(1,1,ncol=(length(occ_mammals_list[[site_sp[1]]])))
  colnames(Sim_commu) <- occ_mammals_list[[site_sp[1]]]
  
  #Compute Ui (Global Uniqueness)
        #Select trait 
  traitsec <- disTraits_birds[,colnames(disTraits_birds) %in% colnames(Sim_commu)]
  traitsec <- traitsec[rownames(traitsec) %in% colnames(Sim_commu),]
  Ui<-uniqueness(Sim_commu,disTraits_birds)
  rownames(Ui)<-Ui[,1]
  #Compute Di (Global distinctiveness)
  Di<-t(distinctiveness(Sim_commu,disTraits_birds))
  colnames(Di)<-"Di"
  
###Birds Mean 
  
  cond <- lapply(l, function(x) length(x) > 3)
  
  
  
  
  
  


LDi_temp <- apply(FR_birds$FR,2,mean,na.rm=TRUE)
LDi <- data.frame(species=names(LDi_temp),LDi=as.numeric(LDi_temp))

GDi <- data.frame(species=FR_birds$FR$species,GDi=FR_birds$FR$Di)

Di__birds <- merge(LDi,GDi)

a <- ggplot(Di__birds, aes(x=GDi,y=LDi))+geom_point() +ggtitle('Birds')

###Birds Var  
LDi_temp <- apply(FR_birds_$Di,2,sd,na.rm=TRUE)
LDi_sd <- data.frame(species=names(LDi_temp),LDi_sd=as.numeric(LDi_temp))

GDi <- data.frame(species=FR_birds_$FR$species,GDi=FR_birds_$FR$Di)

Di__birds <- merge(LDi_sd,GDi)

b <- ggplot(Di__birds, aes(x=GDi,y=LDi_sd))+geom_point() +ggtitle('Birds')


###mammals Mean 
LDi_temp <- apply(FR_mammals_$Di,2,mean,na.rm=TRUE)
LDi <- data.frame(species=names(LDi_temp),LDi=as.numeric(LDi_temp))

GDi <- data.frame(species=FR_mammals_$FR$species,GDi=FR_mammals_$FR$Di)

Di__mammals <- merge(LDi,GDi)

c <- ggplot(Di__mammals, aes(x=GDi,y=LDi))+geom_point() +ggtitle('Mammals')

###mammals Var  
LDi_temp <- apply(FR_mammals_$Di,2,sd,na.rm=TRUE)
LDi_sd <- data.frame(species=names(LDi_temp),LDi_sd=as.numeric(LDi_temp))

GDi <- data.frame(species=FR_mammals_$FR$species,GDi=FR_mammals_$FR$Di)

Di__mammals <- merge(LDi_sd,GDi)

d <- ggplot(Di__mammals, aes(x=GDi,y=LDi_sd))+geom_point() +ggtitle('Mammals')


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








