############
#DATA UPLOAD
#Metadata for the traits : http://www.esapubs.org/archive/ecol/E095/178/metadata.php
############

rm(list=ls(all=TRUE)) 

library(parallel)
library(ggplot2)

source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NM")

#BIRDS Traits ----

    birdstrait_all<- read.csv(file.path(data_dir,"birds","BirdFuncDat.csv"),sep=";")    

    birdstrait_all$Scientific <- gsub(" ","_",birdstrait_all$Scientific)
    birdstrait_all$English <- gsub(" ","_",birdstrait_all$English)
    
    birdstrait_all$Diet.Inv <- as.numeric(birdstrait_all$Diet.Inv)
    birdstrait_all$Diet.Vend <- as.numeric(birdstrait_all$Diet.Vend)
    birdstrait_all$Diet.Vect <- as.numeric(birdstrait_all$Diet.Vect)
    birdstrait_all$Diet.Vfish <- as.numeric(birdstrait_all$Diet.Vfish)
    birdstrait_all$Diet.Vunk <- as.numeric(birdstrait_all$Diet.Vunk)
    birdstrait_all$Diet.Scav <- as.numeric(birdstrait_all$Diet.Scav)
    birdstrait_all$Diet.Fruit <- as.numeric(birdstrait_all$Diet.Fruit)
    birdstrait_all$Diet.Nect <- as.numeric(birdstrait_all$Diet.Nect)
    birdstrait_all$Diet.Seed <- as.numeric(birdstrait_all$Diet.Seed)
    birdstrait_all$Diet.PlantO <- as.numeric(birdstrait_all$Diet.PlantO)
    birdstrait_all$Diet.5Cat <- as.factor(birdstrait_all$Diet.5Cat)
    birdstrait_all$ForStrat.watbelowsurf <- as.numeric(birdstrait_all$ForStrat.watbelowsurf)
    birdstrait_all$ForStrat.wataroundsurf <- as.numeric(birdstrait_all$ForStrat.wataroundsurf)
    birdstrait_all$ForStrat.ground <- as.numeric(birdstrait_all$ForStrat.ground)
    birdstrait_all$ForStrat.understory <- as.numeric(birdstrait_all$ForStrat.understory)
    birdstrait_all$ForStrat.midhigh <- as.numeric(birdstrait_all$ForStrat.midhigh)
    birdstrait_all$ForStrat.canopy <- as.numeric(birdstrait_all$ForStrat.canopy)
    birdstrait_all$ForStrat.aerial <- as.numeric(birdstrait_all$ForStrat.aerial)
    birdstrait_all$PelagicSpecialist <- as.factor(birdstrait_all$PelagicSpecialist)
    birdstrait_all$Nocturnal <- as.factor(birdstrait_all$Nocturnal)
    birdstrait_all$BodyMass.Value <- as.numeric(birdstrait_all$BodyMass.Value)
    
    interest.var.birds <- c("SpecID","Scientific","Diet.Inv","Diet.Vend","Diet.Vect","Diet.Vfish","Diet.Vunk","Diet.Scav","Diet.Fruit","Diet.Nect","Diet.Seed",
                            "Diet.PlantO","Diet.5Cat","ForStrat.watbelowsurf","ForStrat.wataroundsurf","ForStrat.ground","ForStrat.understory",
                            "ForStrat.midhigh","ForStrat.canopy","ForStrat.aerial","PelagicSpecialist","Nocturnal","BodyMass.Value")
    
    birdstrait <- birdstrait_all[,interest.var.birds]
    
    #check uniqueness of SpecID
      duplicated2 <- function(x){
        if (sum(dup <- duplicated(x))==0)
          return(dup)
        if (class(x) %in% c("data.frame","matrix"))
          duplicated(rbind(x[dup,],x))[-(1:sum(dup))]
        else duplicated(c(x[dup],x))[-(1:sum(dup))]
      }
      birdstrait[duplicated2(birdstrait$SpecID),]
      birdstrait <- birdstrait[-c(9994,9995),]
    
    #check if there is NA values 
      lapply(interest.var.birds,function(i)  sum(is.na(birdstrait[,i])))
      
    # keep only the species present in the occ_mamals file 
      load(file=paste0(results_dir,"/birds/occ_birds.RData"))
      birdstrait <- birdstrait[birdstrait$SpecID %in% colnames(occ_birds),]
      
    # Check if the species present in the occ_mamals file have all presence>0 and remove them 
      sumbird <- apply(occ_birds,2,sum)
      sum(sumbird==0)  # ok for all 
      
    #save ID
      birds_id <- birdstrait[,c("SpecID","Scientific")]
      write.csv2(birds_id,file=file.path(results_dir,"birds/Birds_ID.csv"),row.names = FALSE)
    
    #rownames  
      
      rownames(birdstrait) <- paste0("ID_",birdstrait$SpecID)
      birdstrait <- birdstrait[,c(-1,-2)]
    
    #check normality for body mass 
      hist(log(birdstrait$BodyMass.Value))
      
    #test correlation among traits
      test.corel(birdstrait)
      
    #Delete correlated traits
         ##We deleted Diet.5Cat because Rcramer = 0.86 with Diet.Inv
      birdstrait<-birdstrait[,!colnames(birdstrait)%in%"Diet.5Cat"]
      
    #save traits 
      save(birdstrait,file=file.path(results_dir,"birds/birdstrait.RData"))
      
#----
 
#MAMALS Traits ----
    
    mamalstrait_all<- read.csv(file.path(data_dir,"mamals","MamFuncDat.csv"),sep=";") 
    
    mamalstrait_all$Scientific <- gsub(" ","_",mamalstrait_all$Scientific)
    
    mamalstrait_all$Diet.Inv <- as.numeric(mamalstrait_all$Diet.Inv)
    mamalstrait_all$Diet.Vend <- as.numeric(mamalstrait_all$Diet.Vend)
    mamalstrait_all$Diet.Vect <- as.numeric(mamalstrait_all$Diet.Vect)
    mamalstrait_all$Diet.Vfish <- as.numeric(mamalstrait_all$Diet.Vfish)
    mamalstrait_all$Diet.Vunk <- as.numeric(mamalstrait_all$Diet.Vunk)
    mamalstrait_all$Diet.Scav <- as.numeric(mamalstrait_all$Diet.Scav)
    mamalstrait_all$Diet.Fruit <- as.numeric(mamalstrait_all$Diet.Fruit)
    mamalstrait_all$Diet.Nect <- as.numeric(mamalstrait_all$Diet.Nect)
    mamalstrait_all$Diet.Seed <- as.numeric(mamalstrait_all$Diet.Seed)
    mamalstrait_all$Diet.PlantO <- as.numeric(mamalstrait_all$Diet.PlantO)
    mamalstrait_all$ForStrat.Value<- as.factor(mamalstrait_all$ForStrat.Value)
    mamalstrait_all$Activity.Nocturnal<- as.factor(mamalstrait_all$Activity.Nocturnal)
    mamalstrait_all$Activity.Crepuscular<- as.factor(mamalstrait_all$Activity.Crepuscular)
    mamalstrait_all$Activity.Diurnal<- as.factor(mamalstrait_all$Activity.Diurnal)
    mamalstrait_all$BodyMass.Value<- as.numeric(mamalstrait_all$BodyMass.Value)
    
    interest.var.mamals <- c("MSW3_ID","Scientific","Diet.Inv","Diet.Vend","Diet.Vect","Diet.Vfish","Diet.Vunk","Diet.Scav","Diet.Fruit","Diet.Nect","Diet.Seed",
                             "Diet.PlantO","ForStrat.Value","Activity.Nocturnal","Activity.Crepuscular","Activity.Diurnal","BodyMass.Value")
    
    mamalstrait <- mamalstrait_all[,interest.var.mamals]
    
    #check uniqueness of SpecID
      duplicated2 <- function(x){
        if (sum(dup <- duplicated(x))==0)
          return(dup)
        if (class(x) %in% c("data.frame","matrix"))
          duplicated(rbind(x[dup,],x))[-(1:sum(dup))]
        else duplicated(c(x[dup],x))[-(1:sum(dup))]
      }
      mamalstrait[duplicated2(mamalstrait$MSW3_ID),]
    
    #check if there is NA values 
      lapply(interest.var.mamals,function(i)  sum(is.na(mamalstrait[,i])))
      
    # keep only the species present in the occ_mamals file 
      load(file=paste0(results_dir,"/mamals/occ_mamals.RData"))
      mamalstrait <- mamalstrait[mamalstrait$MSW3_ID %in% colnames(occ_mamals),]
      
    # Check if the species present in the occ_mamals file have all presence>0 and remove them 
      summam <- apply(occ_mamals,2,sum)
      sum(summam==0)  # ID_376  ID_437 ID_1243 ID_4066 have occ=0 in all sites 
      occ_mamals <- occ_mamals[,summam>0]
      mamalstrait <- mamalstrait[mamalstrait$MSW3_ID %in% colnames(occ_mamals),]
      
    #save ID
      mamals_id <- mamalstrait[,c("MSW3_ID","Scientific")]
      write.csv2(mamals_id,file=file.path(results_dir,"mamals/Mamals_ID.csv"),row.names = TRUE)
      
    #rownames  
      rownames(mamalstrait) <- paste0("ID_",mamalstrait$MSW3_ID)
      mamalstrait <- mamalstrait[,c(-1,-2)]
    
    #check normality for body mass 
      hist(log(mamalstrait$BodyMass.Value))  
      
    #test correlation among traits
      test.corel(mamalstrait)
 
    #save traits 
      save(mamalstrait,file=file.path(results_dir,"mamals/mamalstrait.RData"))
      
      
#----
