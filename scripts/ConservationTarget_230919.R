rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NL")
library(taxize)

#TODO ANNOTE  SCRIPT
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
load(file=file.path(results_dir,"birds","50km","FR_birds.RData"))
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


save(data_DR_birds, file=file.path(results_dir,"birds/data_DR_birds.RData"))
save(data_DR_mammals, file=file.path(results_dir,"mammals/data_DR_mammals.RData"))


load(file=file.path(results_dir,"birds/data_DR_birds.RData"))
load(file=file.path(results_dir,"mammals/data_DR_mammals.RData"))


### Quantiles of SR and gap analyses

          
       # Mammals -----
        Target_mammals<-data.frame(table(unlist(occ_mammals_list)))
        colnames(Target_mammals) <- c("ID","SR")
        
        Target_mammals$LogSR <- log(Target_mammals$SR) 
        qt=quantile(Target_mammals[,"SR"], probs=c(0.1, 0.9))
        rownames(Target_mammals) <- Target_mammals$ID
        Target_mammals<-Target_mammals[,-1]
        
        MamAllCat<-read.csv2(file=file.path(results_dir,"mammals","50km","MamAllCat.csv"))
        colnames(MamAllCat)[1]<-"Name"
        
        MamAllCat<-merge(MamAllCat,mammalsID,by="Name",all.y=TRUE)
        rownames(MamAllCat)<-MamAllCat$ID
        MamAllCat$PERCENTAGE<-as.numeric(as.character(MamAllCat$PERCENTAGE))
        MamAllCat$PERCENTAGE[is.na(MamAllCat$PERCENTAGE)]<-0
        Target_mammals <- merge (Target_mammals,MamAllCat, by = "row.names" )
        rownames(Target_mammals) <- Target_mammals[,1]
        Target_mammals <- Target_mammals[,c(2,3,5)]
        
        
        Target_mammals[,"TargetExp"] <- target_func(Target_mammals[,"SR"], qt, log=T)
        Target_mammals[,"TargetMet_Percentagecover"] <- 100*(Target_mammals[,"PERCENTAGE"]/Target_mammals[,"TargetExp"])
        Target_mammals <- merge(Target_mammals,data_DR_mammals,by="row.names")
        rownames(Target_mammals) <- Target_mammals[,1]
        Target_mammals <- Target_mammals[,-1]
        #Target_mammals <- na.omit(Target_mammals)

        
        
        ymax=300
        col_br<-c("#00AFBB","#E7B800","orangered")
        Target_mammals_sub <- Target_mammals[((Target_mammals$DR_class=='D25R25') | (Target_mammals$DR_class=='D75R75') | (Target_mammals$DR_class=='AVG')),]
        
        a <- ggplot(Target_mammals_sub, aes(x=DR_class, y=TargetMet_Percentagecover, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
          geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_mammals$TargetMet_Percentagecover,na.rm=T),col="red",linetype="dashed") + 
          labs(x = "DR class",y="Species target achievements")+theme_bw()

        
        # birds -----
        Target_birds<-data.frame(table(unlist(occ_birds_list)))
        colnames(Target_birds) <- c("ID","SR")
        
        Target_birds$LogSR <- log(Target_birds$SR) 
        qt=quantile(Target_birds[,"SR"], probs=c(0.1, 0.9))
        rownames(Target_birds) <- Target_birds$ID
        Target_birds<-Target_birds[,-1]
        
        MamAllCat<-read.csv2(file=file.path(results_dir,"birds","50km","BirdsAllCat.csv"))
        colnames(MamAllCat)[1]<-"Name"
        
        MamAllCat<-merge(MamAllCat,birdsID,by="Name",all.y=TRUE)
        rownames(MamAllCat)<-MamAllCat$ID
        novalue_percentage_birds<-subset(MamAllCat,is.na(MamAllCat$PERCENTAGE))
        save(novalue_percentage_birds, file=file.path(results_dir,"birds/novalue_percentage_birds.RData"))   
          
        MamAllCat$PERCENTAGE<-as.numeric(as.character(MamAllCat$PERCENTAGE))
        MamAllCat$PERCENTAGE[is.na(MamAllCat$PERCENTAGE)]<-0
        Target_birds <- merge (Target_birds,MamAllCat, by = "row.names" )
        rownames(Target_birds) <- Target_birds[,1]
        Target_birds <- Target_birds[,c(2,3,5)]
        
        
        Target_birds[,"TargetExp"] <- target_func(Target_birds[,"SR"], qt, log=T)
        Target_birds[,"TargetMet_Percentagecover"] <- 100*(Target_birds[,"PERCENTAGE"]/Target_birds[,"TargetExp"])
        Target_birds <- merge(Target_birds,data_DR_birds,by="row.names")
        rownames(Target_birds) <- Target_birds[,1]
        Target_birds <- Target_birds[,-1]



ymax=300
col_br<-c("#00AFBB","#E7B800","orangered")
Target_birds_sub <- Target_birds[((Target_birds$DR_class=='D25R25') | (Target_birds$DR_class=='D75R75') | (Target_birds$DR_class=='AVG')),]
b <- ggplot(Target_birds_sub, aes(x=DR_class, y=TargetMet_Percentagecover, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_birds$TargetMet_Percentagecover,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species target achievements")+theme_bw()

grid.arrange(a,b,ncol=2,top = textGrob("Species target achievements" ,gp=gpar(fontsize=20,font=3)))



###Country######
#Analysis of country
####################

load(file=file.path(data_dir,"CountryGrid50km.RData"))
#cells_species <- mclapply(1:nrow(mammalsID),function(i) {grep(mammalsID[i,1], occ_mammals_list)},mc.cores = 3)
#cells_species2 <- mclapply(1:length(cells_species),function(i){ names(occ_mammals_list[cells_species[[i]]])},mc.cores = 3)
#save(cells_species2,file=file.path(results_dir,"mammals","50km","cells_species2.RData"))

country <- data.frame(ID=dataGrid50km$ID, Coundry=dataGrid50km$Country)
country_species <- mclapply(1:length(cells_species2),function(i){unique(country[country$ID %in%  cells_species2[[i]],]$Coundry)},mc.cores = 3)
names(country_species) <-  mammalsID[,1]

country_species_D75R75_mammals <- country_species[names(country_species) %in% rownames(subset(data_DR_mammals,data_DR_mammals$DR_class=="D75R75"))]
country_species_D75R75_mammals<- data.frame(table(unlist(country_species_D75R75_mammals)))
#save(country_species_D75R75_mammals,file=file.path(results_dir,"mammals","50km","country_species_D75R75_mammals.RData"))


cells_species <- mclapply(1:nrow(birdsID),function(i) {grep(birdsID[i,1], occ_birds_list)},mc.cores = 3)
cells_species_birds <- mclapply(1:length(cells_species),function(i){ names(occ_birds_list[cells_species[[i]]])},mc.cores = 3)
save(cells_species_birds,file=file.path(results_dir,"birds","50km","cells_species_birds.RData"))


country <- data.frame(ID=dataGrid50km$ID, country=dataGrid50km$Country)
country_species <- mclapply(1:length(cells_species2),function(i){unique(country[country$ID %in%  cells_species2[[i]],]$country)},mc.cores = 3)
names(country_species) <-  birdsID[,1]

country_species_D75R75_birds <- country_species[names(country_species) %in% rownames(subset(data_DR_birds,data_DR_birds$DR_class=="D75R75"))]
country_species_D75R75_birds<- data.frame(table(unlist(country_species_D75R75_birds)))
#save(country_species_D75R75_birds,file=file.path(results_dir,"birds","50km","country_species_D75R75_birds.RData"))





# load the library
library(forcats)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
#Plot -----
load(file=file.path(results_dir,"birds","50km","country_species_D75R75_birds.RData"))
country_species_D75R75_birds <- subset(country_species_D75R75_birds,country_species_D75R75_birds$Freq >0)
colnames(country_species_D75R75_birds) <- c("country","Nb_SP")
country_species_D75R75_birds <- merge(country_species_D75R75_birds,HDI,by="country")
country_species_D75R75_birds <- country_species_D75R75_birds %>%
  mutate(color = ifelse(country_species_D75R75_birds$Nb_SP >= 10, "red", "black"))

country_species_D75R75_birds <-country_species_D75R75_birds[rev(order(country_species_D75R75_birds$HDI)),]

a <- ggplot(country_species_D75R75_birds,aes(x=HDI, y=reorder(country,HDI),  fill=Nb_SP,colour=Nb_SP,size=Nb_SP)) +
  geom_point(alpha=0.5, shape=21) +
  scale_size(range = c(.1, 10),breaks=c(1,5,10,25,50,75,93), name="Number of rare species", guide="legend") +
  scale_fill_viridis(discrete=FALSE, guide=FALSE, option="C",name="Number of rare species" ) +
  scale_colour_viridis(discrete=FALSE, guide=FALSE, option="C") +
  theme_bw()+theme(legend.position="right",panel.grid.minor = element_blank(), axis.text.y = element_text(size=5.5,color = rev(country_species_D75R75_birds$color)))+ 
  ylab("Country") +
  xlab("HDI") 

load(file=file.path(results_dir,"birds","50km","country_species_D75R75_birds.RData"))
country_species_D75R75_birds <- subset(country_species_D75R75_birds,country_species_D75R75_birds$Freq >0)
colnames(country_species_D75R75_birds) <- c("country","Nb_SP")
country_species_D75R75_birds <- merge(country_species_D75R75_birds,Conflict,by="country")
country_species_D75R75_birds <- na.omit(country_species_D75R75_birds)
country_species_D75R75_birds <- country_species_D75R75_birds %>%
  mutate(color = ifelse(country_species_D75R75_birds$Nb_SP >= 10, "red", "black"))

country_species_D75R75_birds <- country_species_D75R75_birds[rev(order(country_species_D75R75_birds$conflict)),]

b<- ggplot(country_species_D75R75_birds,aes(x=conflict, y=reorder(country,conflict),  fill=Nb_SP,colour=Nb_SP,size=Nb_SP)) +
  geom_point(alpha=0.5, shape=21) +
  scale_size(range = c(.1, 10), breaks=c(1,5,10,25,50,75,93),name="Number of rare species",guide="legend") +
  scale_fill_viridis(discrete=FALSE, guide=FALSE, option="C",name="Number of rare species" ) +
  scale_colour_viridis(discrete=FALSE, guide=FALSE, option="C") +
  theme_bw()+theme(legend.position="right",panel.grid.minor = element_blank(), axis.text.y = element_text(size=5.5,color = rev(country_species_D75R75_birds$color)))+ 
  ylab("Country") +
  xlab("Number of conflicts") 

#Mammals ---
load(file=file.path(results_dir,"mammals","50km","country_species_D75R75_mammals.RData"))
country_species_D75R75_mammals <- subset(country_species_D75R75_mammals,country_species_D75R75_mammals$Freq >0)
colnames(country_species_D75R75_mammals) <- c("country","Nb_SP")
country_species_D75R75_mammals <- merge(country_species_D75R75_mammals,HDI,by="country")
country_species_D75R75_mammals <- country_species_D75R75_mammals %>%
  mutate(color = ifelse(country_species_D75R75_mammals$Nb_SP >= 5, "red", "black"))

country_species_D75R75_mammals <-country_species_D75R75_mammals[rev(order(country_species_D75R75_mammals$HDI)),]

c <- ggplot(country_species_D75R75_mammals,aes(x=HDI, y=reorder(country,HDI),  fill=Nb_SP,colour=Nb_SP,size=Nb_SP)) +
  geom_point(alpha=0.5, shape=21) +
  scale_size(range = c(.1, 10), breaks=c(1,5,10,15,20,30,48), name="Number of rare species",guide="legend") +
  scale_fill_viridis(discrete=FALSE, guide=FALSE, option="C",name="Number of rare species" ) +
  scale_colour_viridis(discrete=FALSE, guide=FALSE, option="C") +
  theme_bw()+theme(legend.position="right",panel.grid.minor = element_blank(), axis.text.y = element_text(size=5.5,color = rev(country_species_D75R75_mammals$color)))+ 
  ylab("Country") +
  xlab("HDI") 

load(file=file.path(results_dir,"mammals","50km","country_species_D75R75_mammals.RData"))
country_species_D75R75_mammals <- subset(country_species_D75R75_mammals,country_species_D75R75_mammals$Freq >0)
colnames(country_species_D75R75_mammals) <- c("country","Nb_SP")
country_species_D75R75_mammals <- merge(country_species_D75R75_mammals,Conflict,by="country")
country_species_D75R75_mammals <- na.omit(country_species_D75R75_mammals)
country_species_D75R75_mammals <- country_species_D75R75_mammals %>%
  mutate(color = ifelse(country_species_D75R75_mammals$Nb_SP >= 5, "red", "black"))

country_species_D75R75_mammals <- country_species_D75R75_mammals[rev(order(country_species_D75R75_mammals$conflict)),]

d<- ggplot(country_species_D75R75_mammals,aes(x=conflict, y=reorder(country,conflict),  fill=Nb_SP,colour=Nb_SP,size=Nb_SP)) +
  geom_point(alpha=0.5, shape=21) +
  scale_size(range = c(.1, 10), breaks=c(1,5,10,15,20,30,48), name="Number of rare species", guide="legend") +
  scale_fill_viridis(discrete=FALSE, guide=FALSE, option="C",name="Number of rare species" ) +
  scale_colour_viridis(discrete=FALSE, guide=FALSE, option="C") +
  theme_bw()+theme(legend.position="right",panel.grid.minor = element_blank(), axis.text.y = element_text(size=5.5,color = rev(country_species_D75R75_mammals$color)))+ 
  ylab("Country") +
  xlab("Number of conflicts") 



#Human foot print
####################

#PER CELLS-----
load(file=file.path(data_dir,"HumanFootprint/dataHF.Rdata"))
load(file=file.path(results_dir,"mammals","50km","funk_mammals.RData"))
load(file=file.path(results_dir,"birds","50km","funk_birds.RData"))

dataHF$ResHF<-as.numeric(as.character(dataHF$ResHF))

Humanfoot_rarety_mammals <- merge(dataHF,funk_mammals,by.x="ID",by.y="cell")
Humanfoot_rarety_birds <- merge(dataHF,funk_birds,by.x="ID",by.y="cell")
Humanfoot_rarety_mammals_AVG <- subset(Humanfoot_rarety_mammals,Humanfoot_rarety_mammals$AVG>0)
Humanfoot_rarety_mammals_D75R75 <- subset(Humanfoot_rarety_mammals,Humanfoot_rarety_mammals$D75R75>0)
Humanfoot_rarety_mammals_D25R25 <- subset(Humanfoot_rarety_mammals,Humanfoot_rarety_mammals$D25R25>0)
Humanfoot_rarety_mammals_NO_D75R75 <- subset(Humanfoot_rarety_mammals,Humanfoot_rarety_mammals$D75R75==0)

mean(Humanfoot_rarety_mammals_D25R25$ResHF,na.rm=T)
sd(Humanfoot_rarety_mammals_D25R25$ResHF,na.rm=T)


HF_mammals <- rbind(data.frame(Value = Humanfoot_rarety_mammals_AVG$ResHF, DR_class = rep("AVG",nrow(Humanfoot_rarety_mammals_AVG)),
                               Threats = rep("Human FootPrint",nrow(Humanfoot_rarety_mammals_AVG)), Taxa =rep("mammals",nrow(Humanfoot_rarety_mammals_AVG))),
                    data.frame(Value = Humanfoot_rarety_mammals_D75R75$ResHF, DR_class = rep("D75R75",nrow(Humanfoot_rarety_mammals_D75R75)),
                               Threats = rep("Human FootPrint",nrow(Humanfoot_rarety_mammals_D75R75)), Taxa =rep("mammals",nrow(Humanfoot_rarety_mammals_D75R75))),
                    data.frame(Value = Humanfoot_rarety_mammals_D25R25$ResHF, DR_class = rep("D25R25",nrow(Humanfoot_rarety_mammals_D25R25)),
                               Threats = rep("Human FootPrint",nrow(Humanfoot_rarety_mammals_D25R25)), Taxa =rep("mammals",nrow(Humanfoot_rarety_mammals_D25R25))))

Humanfoot_rarety_birds_AVG <- subset(Humanfoot_rarety_birds,Humanfoot_rarety_birds$AVG>0)
Humanfoot_rarety_birds_D75R75 <- subset(Humanfoot_rarety_birds,Humanfoot_rarety_birds$D75R75>0)
Humanfoot_rarety_birds_D25R25 <- subset(Humanfoot_rarety_birds,Humanfoot_rarety_birds$D25R25>0)
Humanfoot_rarety_birds_NO_D75R75 <- subset(Humanfoot_rarety_birds,Humanfoot_rarety_birds$D75R75==0)

HF_birds <- rbind(data.frame(Value = Humanfoot_rarety_birds_AVG$ResHF, DR_class = rep("AVG",nrow(Humanfoot_rarety_birds_AVG)),
                               Threats = rep("Human FootPrint",nrow(Humanfoot_rarety_birds_AVG)), Taxa =rep("birds",nrow(Humanfoot_rarety_birds_AVG))),
                    data.frame(Value = Humanfoot_rarety_birds_D75R75$ResHF, DR_class = rep("D75R75",nrow(Humanfoot_rarety_birds_D75R75)),
                               Threats = rep("Human FootPrint",nrow(Humanfoot_rarety_birds_D75R75)), Taxa =rep("birds",nrow(Humanfoot_rarety_birds_D75R75))),
                    data.frame(Value = Humanfoot_rarety_birds_D25R25$ResHF, DR_class = rep("D25R25",nrow(Humanfoot_rarety_birds_D25R25)),
                               Threats = rep("Human FootPrint",nrow(Humanfoot_rarety_birds_D25R25)), Taxa =rep("birds",nrow(Humanfoot_rarety_birds_D25R25))))
HF_mammals_birds <- rbind(HF_mammals,HF_birds)
save(HF_mammals_birds, file=file.path(results_dir,"HF_mammals_birds.RData"))

threats <- rbind(threats,HF_mammals,HF_birds)
par(mfrow=c(1,2))
t.test(subset(Humanfoot_rarety,Humanfoot_rarety$D75R75>0)$ResHF,subset(Humanfoot_rarety,Humanfoot_rarety$D75R75==0)$ResHF,main="mammals")


Humanfoot_rarety <- merge(dataHF,funk_birds,by.x="ID",by.y="cell")
boxplot(subset(Humanfoot_rarety,Humanfoot_rarety$D75R75>0)$ResHF,subset(Humanfoot_rarety,Humanfoot_rarety$D25R25>0)$ResHF,main="birds")   




#PER SPECIEs-----

        #Mammals---
          
          load(file=file.path(results_dir,"mammals","50km","funk_mammals.RData"))
          load(file=file.path(results_dir,"mammals","50km","cells_species_mammals.RData"))
          
          load(file=file.path(data_dir,"HumanFootprint/dataHF.Rdata"))
          dataHF$ResHF<-as.numeric(as.character(dataHF$ResHF))
              # Mammals
          load(file=file.path(data_dir,"CountryGrid50km.RData"))
          
          HDI_species_mammals <- mclapply(1:length(cells_species_mammals),function(i){mean(dataGrid50km[dataGrid50km$ID %in%  cells_species_mammals[[i]],]$HDI2017,na.rm=T)},mc.cores = 3)
          names(HDI_species_mammals) <-  mammalsID[,1]
          HDI_species_mammals<-do.call(rbind,HDI_species_mammals)
          
          ConflictCY_species_mammals <- mclapply(1:length(cells_species_mammals),function(i){mean(dataGrid50km[dataGrid50km$ID %in%  cells_species_mammals[[i]],]$ConflictCY,na.rm=T)},mc.cores = 3)
          names(ConflictCY_species_mammals) <-  mammalsID[,1]
          ConflictCY_species_mammals<-do.call(rbind,ConflictCY_species_mammals)
          
          HF_species_mammals <- mclapply(1:length(cells_species_mammals),function(i){mean(dataHF[dataHF$ID %in%  cells_species_mammals[[i]],]$ResHF,na.rm=T)},mc.cores = 3)
          names(HF_species_mammals) <-  mammalsID[,1]
          HF_species_mammals<-do.call(rbind,HF_species_mammals)
          
          #Function to merge
          MyMerge <- function(x, y){
            df    <- merge(x, y, by= "row.names", all.x= T, all.y= F)
            rownames(df)  <- df$Row.names
            df$Row.names  <- NULL
            return(df)
          }
          
          threats_species_mammals <- Reduce(MyMerge, list(HDI_species_mammals, ConflictCY_species_mammals, HF_species_mammals, Target_mammals))
          colnames(threats_species_mammals) <- c("HDI","Conflict","HumanFootPrint","SR", "LogSR", "PERCENTAGE", "TargetExp", "TargetMet_Percentagecover")
          
          threats_species_mammals <- merge(threats_species_mammals,data_DR_mammals,by="row.names")
          rownames(threats_species_mammals) <- threats_species_mammals[,1]
          threats_species_mammals <- threats_species_mammals[,-1]
          threats_species_mammals <- merge(threats_species_mammals,mammalsID,by.x="row.names",by.y="ID")
          rownames(threats_species_mammals) <- threats_species_mammals[,1]
          threats_species_mammals <- threats_species_mammals[,-1]
          
      
          #birds---
          
          load(file=file.path(results_dir,"birds","50km","funk_birds.RData"))
          load(file=file.path(results_dir,"birds","50km","cells_species_birds.RData"))
          load(file=file.path(results_dir,"birds","50km","data_DR_birds.RData"))
          
          rownames(Target_birds) <- Target_birds[,1]
          Target_birds <- Target_birds[,-1]
          
          HDI_species_birds <- mclapply(1:length(cells_species_birds),function(i){mean(dataGrid50km[dataGrid50km$ID %in%  cells_species_birds[[i]],]$HDI2017,na.rm=T)},mc.cores = 3)
          names(HDI_species_birds) <-  birdsID[,1]
          HDI_species_birds<-do.call(rbind,HDI_species_birds)
          
          ConflictCY_species_birds <- mclapply(1:length(cells_species_birds),function(i){mean(dataGrid50km[dataGrid50km$ID %in%  cells_species_birds[[i]],]$ConflictCY,na.rm=T)},mc.cores = 3)
          names(ConflictCY_species_birds) <-  birdsID[,1]
          ConflictCY_species_birds<-do.call(rbind,ConflictCY_species_birds)
          
          HF_species_birds <- mclapply(1:length(cells_species_birds),function(i){mean(dataHF[dataHF$ID %in%  cells_species_birds[[i]],]$ResHF,na.rm=T)},mc.cores = 3)
          names(HF_species_birds) <-  birdsID[,1]
          HF_species_birds<-do.call(rbind,HF_species_birds)
          
          #Function to merge
          MyMerge <- function(x, y){
            df    <- merge(x, y, by= "row.names", all.x= T, all.y= F)
            rownames(df)  <- df$Row.names
            df$Row.names  <- NULL
            return(df)
          }
          
          threats_species_birds <- Reduce(MyMerge, list(HDI_species_birds, ConflictCY_species_birds, HF_species_birds, Target_birds))
          colnames(threats_species_birds) <- c("HDI","Conflict","HumanFootPrint","SR", "LogSR", "PERCENTAGE", "TargetExp", "TargetMet_Percentagecover")
          
          threats_species_birds <- merge(threats_species_birds,data_DR_birds,by="row.names",y.all=T)
          rownames(threats_species_birds) <- threats_species_birds[,1]
          threats_species_birds <- threats_species_birds[,-1]
          threats_species_birds <- merge(threats_species_birds,birdsID,by.x="row.names",by.y="ID")
          rownames(threats_species_birds) <- threats_species_birds[,1]
          threats_species_birds <- threats_species_birds[,-1]
          

threats_per_species<-data.frame(value = c(threats_species_birds$HDI,threats_species_mammals$HDI,
                                        threats_species_birds$Conflict,threats_species_mammals$Conflict,
                                        threats_species_birds$HumanFootPrint,threats_species_mammals$HumanFootPrint,
                                        threats_species_birds$TargetMet_Percentagecover,threats_species_mammals$TargetMet_Percentagecover),
                                
                                DR_class = rep(c(as.character(threats_species_birds$DR_class),as.character(threats_species_mammals$DR_class)),4),
                                
                                Threats = c(rep("meanHDI",length(c(threats_species_birds$HDI,threats_species_mammals$HDI))),
                                            rep("meanConflict",length(c(threats_species_birds$Conflict,threats_species_mammals$Conflict))),
                                            rep("Human FootPrint",length(c(threats_species_birds$HumanFootPrint,threats_species_mammals$HumanFootPrint))),
                                            rep("Protection Target",length(c(threats_species_birds$TargetMet_Percentagecover,threats_species_mammals$TargetMet_Percentagecover)))),
                                
                                Taxa = rep(c(rep("birds",nrow(threats_species_birds)),rep("mammals",nrow(threats_species_mammals))),4))


save(threats_per_species,file=file.path(results_dir,"threats_per_species"))


load("~/Documents/Postdoc MARBEC/FREE/RALLL_R/FUNCRARITY/data/threats_cc_with_id.RData")
climate_birds <- subset(threats_cc, threats_cc$Taxa=="birds")
climate_birds_2040 <- subset(threats_cc, threats_cc$Horizon=="2041-2060")
climate_birds_2060 <- subset(threats_cc, threats_cc$Horizon=="2061-2080")

climate_birds_2040 <- climate_birds_2040[,c(1,2)]
climate_birds_2060 <- climate_birds_2060[,c(1,2)]
threats_species_birds_cc <-merge (threats_species_birds,climate_birds_2040,by.x="Row.names",by.y="ID",all.x=T)
threats_species_birds_cc <-merge (threats_species_birds_cc,climate_birds_2060,by.x="Row.names",by.y="ID",all.x=T)
colnames(threats_species_birds_cc)[14] <- "2041-2060"
colnames(threats_species_birds_cc)[15] <- "2061-2080"


climate_mammals <- subset(threats_cc, threats_cc$Taxa=="mammals")
climate_mammals_2040 <- subset(threats_cc, threats_cc$Horizon=="2041-2060")
climate_mammals_2060 <- subset(threats_cc, threats_cc$Horizon=="2061-2080")

climate_mammals_2040 <- climate_mammals_2040[,c(1,2)]
climate_mammals_2060 <- climate_mammals_2060[,c(1,2)]
threats_species_mammals_cc <-merge (threats_species_mammals,climate_mammals_2040,by.x="row.names",by.y="ID",all.x=T)
threats_species_mammals_cc <-merge (threats_species_mammals_cc,climate_mammals_2060,by.x="Row.names",by.y="ID",all.x=T)
colnames(threats_species_mammals_cc)[16] <- "2041-2060"
colnames(threats_species_mammals_cc)[17] <- "2061-2080"


## NEED TO COMPUTE CENTROID PTO HAVE THE MOST IMPACTED
rare_mammals <- subset(threats_species_mammals_cc,threats_species_mammals_cc$DR_class=="D75R75")
rownames(rare_mammals) <- rare_mammals[,1]

rare_mammals_acp <- dudi.pca(na.omit(threats_species_mammals_cc[,c(2:4,9,16,17)]))
fviz_pca_ind(rare_mammals_acp, axes = c(1, 2), geom = c("point", "text"),
             label = "all", invisible = "none", labelsize = 4,
             pointsize = 2, habillage = "none",
             addEllipses = FALSE, ellipse.level = 0.95, 
             col.ind = "black", col.ind.sup = "blue", alpha.ind = 1,
             select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
             jitter = list(what = "label", width = NULL, height = NULL))
# Graphique des variables
fviz_pca_var(rare_mammals_acp, axes = c(1, 2), geom = c("arrow", "text"),
             label = "all", invisible = "none", labelsize = 4,
             col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
             col.circle = "grey70",
             select.var = list(name =NULL, cos2 = NULL, contrib = NULL),
             jitter = list(what = "label", width = NULL, height = NULL))
# Biplot des individus et des variables
fviz_pca_biplot(rare_mammals_acp, axes = c(1, 2), geom = c("point", "text"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 2,
                habillage = "none", addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 1,
                col.var = "steelblue", alpha.var = 1, col.quanti.sup = "blue",
                col.circle = "grey70", 
                select.var = list(name = NULL, cos2 = NULL, contrib= NULL), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL))


rare_mammals_less_prot<-rare_mammals[order(rare_mammals[,"TargetMet_Percentagecover"],decreasing=F)[1:50],]
rare_mammals_more_HF<-rare_mammals[order(rare_mammals[,"HumanFootPrint"],decreasing=T)[1:50],]
rare_mammals_less_HDI<-rare_mammals[order(rare_mammals[,"HDI"],decreasing=F)[1:50],]
rare_mammals_more_conflict<-rare_mammals[order(rare_mammals[,"Conflict"],decreasing=T)[1:50],]
rare_mammals_climate<-  rare_mammals[order(rare_mammals[,"2041-2060"],decreasing=F)[1:50],]

list_of_data = list(rare_mammals_more_HF,rare_mammals_less_prot,rare_mammals_less_HDI,rare_mammals_more_conflict,rare_mammals_climate)
common_names = Reduce(intersect, lapply(list_of_data, row.names))

list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,]})
list_of_data







rare_birds <- subset(threats_species_birds,threats_species_birds$DR_class=="D75R75")


rare_birds_acp <- dudi.pca(na.omit(rare_birds_threats[,c(1:3,8)]))
fviz_pca_ind(rare_birds_acp, axes = c(1, 2), geom = c("point", "text"),
             label = "all", invisible = "none", labelsize = 4,
             pointsize = 2, habillage = "none",
             addEllipses = FALSE, ellipse.level = 0.95, 
             col.ind = "black", col.ind.sup = "blue", alpha.ind = 1,
             select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
             jitter = list(what = "label", width = NULL, height = NULL))
# Graphique des variables
fviz_pca_var(rare_birds_acp, axes = c(1, 2), geom = c("arrow", "text"),
             label = "all", invisible = "none", labelsize = 4,
             col.var = "black", alpha.var = 1, col.quanti.sup = "blue",
             col.circle = "grey70",
             select.var = list(name =NULL, cos2 = NULL, contrib = NULL),
             jitter = list(what = "label", width = NULL, height = NULL))
# Biplot des individus et des variables
fviz_pca_biplot(rare_birds_acp, axes = c(1, 2), geom = c("point", "text"),
                label = "all", invisible = "none", labelsize = 4, pointsize = 2,
                habillage = "none", addEllipses = FALSE, ellipse.level = 0.95,
                col.ind = "black", col.ind.sup = "blue", alpha.ind = 1,
                col.var = "steelblue", alpha.var = 1, col.quanti.sup = "blue",
                col.circle = "grey70", 
                select.var = list(name = NULL, cos2 = NULL, contrib= NULL), 
                select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                jitter = list(what = "label", width = NULL, height = NULL))


rare_birds_less_prot<-rare_birds[order(rare_birds[,"TargetMet_Percentagecover"],decreasing=F)[1:50],]
rare_birds_more_HF<-rare_birds[order(rare_birds[,"HumanFootPrint"],decreasing=T)[1:50],]
rare_birds_less_HDI<-rare_birds[order(rare_birds[,"HDI"],decreasing=F)[1:50],]
rare_birds_more_conflict<-rare_birds[order(rare_birds[,"Conflict"],decreasing=T)[1:50],]
rare_birds_climate<-  climate[(climate$DR_class=="D75R75") | (climate$Taxa=="birds"),]
rare_birds_climate<-  rare_birds_climate[order(rare_birds_climate[,"Value"],decreasing=F)[1:50],]


test <- subset(rare_birds_acp$li,rare_birds_acp$li$Axis1<0.1 & rare_birds_acp$li$Axis1>(-0.1))
test2 <- subset(test,test$Axis2<0.1 & test$Axis2>(-0.1))
test3 <-  threats_species_birds
test3[test3[,1] %in% rownames(test2),]


list_of_data = list(rare_birds_more_HF,rare_birds_less_prot)
common_names = Reduce(intersect, lapply(list_of_data, row.names))

list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,]})
list_of_data

threats_cc_with_id.RData.                                                                                                                                                                                                                  -0.005417, 0.000522, -0.003114), PCA.7 = c(-0.056734, -0.007418, 
                                                                                                                                                                                                                                                                                                                                                                      -0.01043, -0.006961, -0.006006), PCA.8 = c(0.005189, 0.008031, 
                                                                                                                                                                                                                                                                                                                                                                                                                 -0.002979, 0.000743, 0.006276), PCA.9 = c(0.008169, -0.000265, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                           0.010893, 0.003233, 0.007316)), .Names = c("Sample", "PCA.1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "PCA.2", "PCA.3", "PCA.4", "PCA.5", "PCA.6", "PCA.7", "PCA.8", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "PCA.9"), row.names = c(NA, 5L), class = "data.frame")









#HDI
####################
Conflict<- unique(data.frame(conflict=dataGrid50km$ConflictCY, country=dataGrid50km$Country))

HDI <- unique(data.frame(HDI=dataGrid50km$HDI2017, country=dataGrid50km$Country))
country <- data.frame(ID=dataGrid50km$ID, Coundry=dataGrid50km$Country)
country_rarety <- cbind(country,funk_birds$D75R75,funk_mammals$D75R75)
colnames(country_rarety) <- c("ID","country","birds","mammals")

country_rarety$birds0_1 <- funk_birds$D75R75
country_rarety$mammas0_1 <- funk_mammals$D75R75

country_rarety$birds0_1[country_rarety$birds0_1>0] <-1
country_rarety$mammals0_1[country_rarety$mammas0_1>0] <-1

country_important_birds <- data.frame(xtabs(birds0_1~country,data=country_rarety))
country_important_birds <- subset(country_important_birds,country_important_birds$Freq>0)
boxplot(test$HDI,test2$HDI)

country_important_mammals <- data.frame(xtabs(mammals0_1~country,data=country_rarety))
country_important_mammals <- subset(country_important_mammals,country_important_mammals$Freq>0)

test <- merge(country_important_mammals,HDI,by="country")
test2 <- HDI[!HDI$country %in% country_important_mammals$country,]

boxplot(test$HDI,test2$HDI)
t.test(test$HDI,test2$HDI)
test <- merge(country_important_birds,HDI,by="country")
test2 <- HDI[!HDI$country %in% country_important_birds$country,]

t.test(test$HDI,test2$HDI)




















#####################OLD STUFF

load(file=file.path(results_dir,"mammals","50km","mammals_PA.RData"))


Target_mammals <- merge (Target_mammals,mammals_PA, by = "row.names" )
rownames(Target_mammals) <- Target_mammals[,1]
Target_mammals <- Target_mammals[,c(2,3,5,6,15)]


Target_mammals[,"TargetExp"] <- target_func(Target_mammals[,"SR"], qt, log=T)
Target_mammals[,"TargetMet_Percentagecover"] <- 100*(Target_mammals[,"Percentagecover"]/Target_mammals[,"TargetExp"])
Target_mammals <- merge(Target_mammals,data_DR_mammals,by="row.names")
Target_mammals <- na.omit(Target_mammals)

load(file=file.path(results_dir,"birds","50km","birds_PA.RData"))
Target_birds <- merge (Target_birds,birds_PA, by = "row.names" )
rownames(Target_birds) <- Target_birds[,1]
Target_birds <- Target_birds[,c(2,3,5,6,15)]


Target_birds[,"TargetExp"] <- target_func(Target_birds[,"SR"], qt, log=T)
Target_birds[,"TargetMet_Percentagecover"] <- 100*(Target_birds[,"Percentagecover"]/Target_birds[,"TargetExp"])
Target_birds <- merge(Target_birds,data_DR_birds,by="row.names")
Target_birds <- na.omit(Target_birds)


#Species target achievement was then measured at the ratio between the percentage of
#range actually covered by the PA system and the target.  
ymax=300
col_br<-rev(heat.colors(5))
a <- ggplot(Target_mammals, aes(x=DR_class, y=TargetMet_Percentagecover, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_mammals$TargetMet_Percentagecover,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species’ target achievements")+theme_bw()

Target_mammals_sub <- Target_mammals[((Target_mammals$DR_class=='D25R25') | (Target_mammals$DR_class=='D75R75') | (Target_mammals$DR_class=='AVG')),]
b <- ggplot(Target_mammals_sub, aes(TargetMet_Percentagecover,fill=DR_class,color=DR_class)) + geom_density(adjust = 1.5,alpha = 0.1) + xlim(0, ymax)+ 
  scale_fill_manual(values=c(col_br[1],col_br[2],col_br[5]))+ scale_color_manual(values=c(col_br[1],col_br[2],"red"))+
  theme(legend.position = c(0.9, 0.8)) + geom_vline(xintercept=mean(Target_mammals_sub$TargetMet_Percentagecover,na.rm=T),col="grey35",linetype="dashed")+
  labs(x = "Species target achievements")+theme_bw()

c <- ggplot(Target_birds, aes(x=DR_class, y=TargetMet_Percentagecover, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_birds$TargetMet_Percentagecover,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species target achievements")+theme_bw()

Target_birds_sub <- Target_birds[((Target_birds$DR_class=='D25R25') | (Target_birds$DR_class=='D75R75') | (Target_birds$DR_class=='AVG')),]
d <- ggplot(Target_birds_sub, aes(TargetMet_Percentagecover,fill=DR_class,color=DR_class)) + geom_density(adjust = 1.5,alpha = 0.1) + xlim(0, ymax)+ 
  scale_fill_manual(values=c(col_br[1],col_br[2],col_br[5]))+ scale_color_manual(values=c(col_br[1],col_br[2],"red"))+
  theme(legend.position = c(0.9, 0.8)) + geom_vline(xintercept=mean(Target_birds_sub$TargetMet_Percentagecover,na.rm=T),col="grey35",linetype="dashed")+
  labs(x = "Species target achievements")+theme_bw()

grid.arrange(b,d,ncol=2,top = textGrob("Species target achievements" ,gp=gpar(fontsize=20,font=3)))

#PLutot que la moyenne faire la somme de la surface protégé





#NEw DATA # NEED TO DIVISE BY 100 THE AREA
load(file=file.path(data_dir,"mammals","50km","CatIaMammals.RData"))
load(file=file.path(data_dir,"mammals","50km","CatIbMammals.RData"))
load(file=file.path(data_dir,"mammals","50km","CatIIMammals.RData"))

CatIaMammals <- CatIaMammals[CatIaMammals$BINOMIAL %in% mammalsID$Name,]
CatIaMammals <- merge(CatIaMammals,mammalsID, by.x="BINOMIAL", by.y= "Name")

CatIbMammals <- CatIbMammals[CatIbMammals$BINOMIAL %in% mammalsID$Name,]
CatIbMammals <- merge(CatIbMammals,mammalsID, by.x="BINOMIAL", by.y= "Name")

CatIIMammals <- CatIIMammals[CatIIMammals$BINOMIAL %in% mammalsID$Name,]
CatIIMammals <- merge(CatIIMammals,mammalsID, by.x="BINOMIAL", by.y= "Name")


CatMammals_all <-  data.frame(BINOMIAL=CatIaMammals$BINOMIAL, TOTAL_AREA=CatIaMammals$TOTAL_AREA, AREA_PRO=apply(data.frame(CatIaMammals$AREA/100,CatIbMammals$AREA/100,CatIIMammals$AREA/100),1,sum))
CatMammals_all$PERCENTAGE <- (CatMammals_all$AREA_PRO/CatMammals_all$TOTAL_AREA)*100
CatMammals_all<- merge(CatMammals_all,mammalsID, by.x="BINOMIAL", by.y= "Name")

Target_mammals_all<- data.frame(ID=CatMammals_all$ID, SR=CatMammals_all$TOTAL_AREA,Percentagecover=CatMammals_all$PERCENTAGE)
Target_mammalsIa <- data.frame(ID=CatIaMammals$ID, SR=CatIaMammals$TOTAL_AREA,Percentagecover=CatIaMammals$PERCENTAGE)
Target_mammalsIb <- data.frame(ID=CatIbMammals$ID, SR=CatIbMammals$TOTAL_AREA,Percentagecover=CatIbMammals$PERCENTAGE)
Target_mammalsII <- data.frame(ID=CatIIMammals$ID, SR=CatIIMammals$TOTAL_AREA,Percentagecover=CatIIMammals$PERCENTAGE)

Target_mammals <- cbind(rbind(Target_mammalsIa,Target_mammalsIb,Target_mammalsII), c(rep("Ia",nrow(Target_mammalsIa)),rep("Ib",nrow(Target_mammalsIb)),rep("II",nrow(Target_mammalsII))))
colnames(Target_mammals)[4] <- "IUCN_CAT"
Target_mammals$LogSR <- log(Target_mammals$SR) 
qt=quantile(Target_mammals[,"SR"], probs=c(0.1, 0.9))
#rownames(Target_mammals) <- Target_mammals$ID
#Target_mammals<-Target_mammals[,-1]
Target_mammals[,"TargetExp"] <- target_func(Target_mammals[,"SR"], qt, log=T)
Target_mammals[,"TargetAchiev"] <- 100*(Target_mammals[,"Percentagecover"]/Target_mammals[,"TargetExp"])
Target_mammals <- merge(Target_mammals,data_DR_mammals,by.x="ID",by.y="row.names")
Target_mammals <- na.omit(Target_mammals)

Target_mammals_all$LogSR <- log(Target_mammals_all$SR) 
qt=quantile(Target_mammals_all[,"SR"], probs=c(0.1, 0.9))
Target_mammals_all[,"TargetExp"] <- target_func(Target_mammals_all[,"SR"], qt, log=T)
Target_mammals_all[,"TargetAchiev"] <- 100*(Target_mammals_all[,"Percentagecover"]/Target_mammals_all[,"TargetExp"])
Target_mammals_all <- merge(Target_mammals_all,data_DR_mammals,by.x="ID",by.y="row.names")
Target_mammals_all <- na.omit(Target_mammals_all)

ymax=100
col_br<-rev(heat.colors(3))
Target_mammals_sub <- Target_mammals[((Target_mammals$DR_class=='D25R25') | (Target_mammals$DR_class=='D75R75') | (Target_mammals$DR_class=='AVG')),]
a <- ggplot(Target_mammals_sub, aes(x=DR_class, y=TargetAchiev, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_mammals_sub$TargetAchiev,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species’ target achievements")+theme_bw()+
facet_wrap("IUCN_CAT")


Target_mammals_all_sub <- Target_mammals_all[((Target_mammals_all$DR_class=='D25R25') | (Target_mammals_all$DR_class=='D75R75') | (Target_mammals_all$DR_class=='AVG')),]
b <- ggplot(Target_mammals_all_sub, aes(x=DR_class, y=TargetAchiev, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_mammals_all_sub$TargetAchiev,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species’ target achievements")+theme_bw()

Target_mammals_all_sub$diff <- Target_mammals_all_sub$TargetAchiev - Target_mammals_all_sub$TargetExp
c <- ggplot(Target_mammals_all_sub, aes(x=DR_class, y=diff, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(-100, 100)) + geom_hline(yintercept=mean(Target_mammals_all_sub$diff,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species’ target achievements")+theme_bw()
#--------
load(file=file.path(data_dir,"birds","50km","CatIaBirds.RData"))
load(file=file.path(data_dir,"birds","50km","CatIbBirds.RData"))
load(file=file.path(data_dir,"birds","50km","CatIIBirds.RData"))

CatIaBirds <- CatIaBirds[CatIaBirds$SCINAME %in% birdsID$Name,]
CatIaBirds <- merge(CatIaBirds,birdsID, by.x="SCINAME", by.y= "Name")

CatIbBirds <- CatIbBirds[CatIbBirds$SCINAME %in% birdsID$Name,]
CatIbBirds <- merge(CatIbBirds,birdsID, by.x="SCINAME", by.y= "Name")

CatIIBirds <- CatIIBirds[CatIIBirds$SCINAME %in% birdsID$Name,]
CatIIBirds <- merge(CatIIBirds,birdsID, by.x="SCINAME", by.y= "Name")


CatBirds_all <-  data.frame(SCINAME=CatIaBirds$SCINAME, TOTAL_AREA=CatIaBirds$TOTAL_AREA, AREA_PRO=apply(data.frame(CatIaBirds$AREA/100,CatIbBirds$AREA/100,CatIIBirds$AREA/100),1,sum))
CatBirds_all$PERCENTAGE <- (CatBirds_all$AREA_PRO/CatBirds_all$TOTAL_AREA)*100
CatBirds_all<- merge(CatBirds_all,birdsID, by.x="SCINAME", by.y= "Name")

Target_birds_all<- data.frame(ID=CatBirds_all$ID, SR=CatBirds_all$TOTAL_AREA,Percentagecover=CatBirds_all$PERCENTAGE)
Target_birdsIa <- data.frame(ID=CatIaBirds$ID, SR=CatIaBirds$TOTAL_AREA,Percentagecover=CatIaBirds$PERCENTAGE)
Target_birdsIb <- data.frame(ID=CatIbBirds$ID, SR=CatIbBirds$TOTAL_AREA,Percentagecover=CatIbBirds$PERCENTAGE)
Target_birdsII <- data.frame(ID=CatIIBirds$ID, SR=CatIIBirds$TOTAL_AREA,Percentagecover=CatIIBirds$PERCENTAGE)

Target_birds <- cbind(rbind(Target_birdsIa,Target_birdsIb,Target_birdsII), c(rep("Ia",nrow(Target_birdsIa)),rep("Ib",nrow(Target_birdsIb)),rep("II",nrow(Target_birdsII))))
colnames(Target_birds)[4] <- "IUCN_CAT"
Target_birds$LogSR <- log(Target_birds$SR) 
qt=quantile(Target_birds[,"SR"], probs=c(0.1, 0.9))
#rownames(Target_birds) <- Target_birds$ID
#Target_birds<-Target_birds[,-1]
Target_birds[,"TargetExp"] <- target_func(Target_birds[,"SR"], qt, log=T)
Target_birds[,"TargetAchiev"] <- 100*(Target_birds[,"Percentagecover"]/Target_birds[,"TargetExp"])
Target_birds <- merge(Target_birds,data_DR_birds,by.x="ID",by.y="row.names")
Target_birds <- na.omit(Target_birds)



Target_birds_all$LogSR <- log(Target_birds_all$SR) 
qt=quantile(Target_birds_all[,"SR"], probs=c(0.1, 0.9))
Target_birds_all[,"TargetExp"] <- target_func(Target_birds_all[,"SR"], qt, log=T)
Target_birds_all[,"TargetAchiev"] <- 100*(Target_birds_all[,"Percentagecover"]/Target_birds_all[,"TargetExp"])
Target_birds_all <- merge(Target_birds_all,data_DR_birds,by.x="ID",by.y="row.names")
Target_birds_all <- na.omit(Target_birds_all)


ymax=200
col_br<-rev(heat.colors(3))
Target_birds_sub <- Target_birds[((Target_birds$DR_class=='D25R25') | (Target_birds$DR_class=='D75R75') | (Target_birds$DR_class=='AVG')),]
a <- ggplot(Target_birds_sub, aes(x=DR_class, y=TargetAchiev, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_birds_sub$TargetAchiev,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species’ target achievements")+theme_bw()+
  facet_wrap("IUCN_CAT")


Target_birds_all_sub <- Target_birds_all[((Target_birds_all$DR_class=='D25R25') | (Target_birds_all$DR_class=='D75R75') | (Target_birds_all$DR_class=='AVG')),]

b <- ggplot(Target_birds_all_sub, aes(x=DR_class, y=TargetAchiev, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_birds_all_sub$TargetAchiev,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species’ target achievements")+theme_bw()
Target_birds_all_sub$diff <- Target_birds_all_sub$TargetAchiev - Target_birds_all_sub$TargetExp
c <- ggplot(Target_birds_all_sub, aes(x=DR_class, y=diff, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(-100, 100)) + geom_hline(yintercept=mean(Target_birds_all_sub$diff,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species’ target achievements")+theme_bw()


ggplot(Target_birds_all_sub, aes(x = LogSR, y = Percentagecover, color = TargetAchiev, size = TargetAchiev)) + 
  scale_color_gradient(low="red", high="green")+
  geom_point(alpha = 0.5) +  facet_wrap("DR_class")+theme_bw()

ggplot(Target_birds_all_sub, aes(x = LogSR, y = TargetAchiev, color =TargetAchiev , size = TargetAchiev)) + 
  scale_color_gradient(low="red", high="green")+
  geom_point(alpha = 0.5) +  facet_wrap("DR_class")+theme_bw()


ggplot(Target_birds_all_sub, aes(x = LogSR, y = TargetAchiev, color =TargetAchiev , size = TargetAchiev)) + 
  scale_color_gradient(low="red", high="green")


library(colorRamps)     # for matlab.like(...)
ggp + scale_color_gradientn(colours=matlab.like(10))

# MAMMALS
Target_mammals<-data.frame(table(unlist(occ_mammals_list)))
colnames(Target_mammals) <- c("ID","SR")

Target_mammals$LogSR <- log(Target_mammals$SR) 
qt=quantile(Target_mammals[,"SR"], probs=c(0.1, 0.9))
rownames(Target_mammals) <- Target_mammals$ID
Target_mammals<-Target_mammals[,-1]

MamAllCat<-read.csv2(file=file.path(results_dir,"mammals","50km","MamAllCat.csv"))
MamAllCat<-merge(MamAllCat,mammalsID,by.x="X...SPECIES", by.y="Name")
rownames(MamAllCat)<-MamAllCat$ID
MamAllCat$PERCENTAGE<-as.numeric(as.character(MamAllCat$PERCENTAGE))

Target_mammals <- merge (Target_mammals,MamAllCat, by = "row.names" )
rownames(Target_mammals) <- Target_mammals[,1]
Target_mammals <- Target_mammals[,c(2,3,5)]


Target_mammals[,"TargetExp"] <- target_func(Target_mammals[,"SR"], qt, log=T)
Target_mammals[,"TargetMet_Percentagecover"] <- 100*(Target_mammals[,"PERCENTAGE"]/Target_mammals[,"TargetExp"])
Target_mammals <- merge(Target_mammals,data_DR_mammals,by="row.names")
Target_mammals <- na.omit(Target_mammals)





Target_mammals_sub <- Target_mammals[((Target_mammals$DR_class=='D25R25') | (Target_mammals$DR_class=='D75R75') | (Target_mammals$DR_class=='AVG')),]
ymax=300
col_br<-c("#00AFBB","#E7B800","orangered")
a <- ggplot(Target_mammals_sub, aes(x=DR_class, y=TargetMet_Percentagecover, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_mammals$TargetMet_Percentagecover,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species target achievements")+theme_bw()



# BIRDS
Target_birds<-data.frame(table(unlist(occ_birds_list)))
colnames(Target_birds) <- c("ID","SR")

Target_birds$LogSR <- log(Target_birds$SR) 
qt=quantile(Target_birds[,"SR"], probs=c(0.1, 0.9))
rownames(Target_birds) <- Target_birds$ID
Target_birds<-Target_birds[,-1]

BirdsAllCat<-read.csv2(file=file.path(results_dir,"birds","50km","BirdsAllCat.csv"))
BirdsAllCat<-merge(BirdsAllCat,birdsID,by.x="X...SPECIES", by.y="Name")
rownames(BirdsAllCat)<-BirdsAllCat$ID
BirdsAllCat$PERCENTAGE<-as.numeric(as.character(BirdsAllCat$PERCENTAGE))

Target_birds <- merge (Target_birds,BirdsAllCat, by = "row.names" )
rownames(Target_birds) <- Target_birds[,1]
Target_birds <- Target_birds[,c(2,3,5)]


Target_birds[,"TargetExp"] <- target_func(Target_birds[,"SR"], qt, log=T)
Target_birds[,"TargetMet_Percentagecover"] <- 100*(Target_birds[,"PERCENTAGE"]/Target_birds[,"TargetExp"])
Target_birds <- merge(Target_birds,data_DR_birds,by="row.names")
Target_birds <- na.omit(Target_birds)


ymax=300
col_br<-c("#00AFBB","#E7B800","orangered")
Target_birds_sub <- Target_birds[((Target_birds$DR_class=='D25R25') | (Target_birds$DR_class=='D75R75') | (Target_birds$DR_class=='AVG')),]
b <- ggplot(Target_birds_sub, aes(x=DR_class, y=TargetMet_Percentagecover, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
  geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(Target_birds$TargetMet_Percentagecover,na.rm=T),col="red",linetype="dashed") + 
  labs(x = "DR class",y="Species target achievements")+theme_bw()


median(subset(Target_birds_sub,Target_birds_sub$DR_class=="D75R75")$TargetMet_Percentagecover)
median(subset(Target_birds_sub,Target_birds_sub$DR_class=="D25R25")$TargetMet_Percentagecover)

median(subset(Target_mammals_sub,Target_mammals_sub$DR_class=="D75R75")$TargetMet_Percentagecover)
median(subset(Target_mammals_sub,Target_mammals_sub$DR_class=="D25R25")$TargetMet_Percentagecover)

grid.arrange(a,b,ncol=2,top = textGrob("Species target achievements" ,gp=gpar(fontsize=20,font=3)))

