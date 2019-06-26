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

          # MAMMALS
          Target_mammals<-data.frame(table(unlist(occ_mammals_list)))
          colnames(Target_mammals) <- c("ID","SR")
          
          Target_mammals$LogSR <- log(Target_mammals$SR) 
          qt=quantile(Target_mammals[,"SR"], probs=c(0.1, 0.9))
          rownames(Target_mammals) <- Target_mammals$ID
          Target_mammals<-Target_mammals[,-1]
          
          # BIRDS
          Target_birds<-data.frame(table(unlist(occ_birds_list)))
          colnames(Target_birds) <- c("ID","SR")
          
          Target_birds$LogSR <- log(Target_birds$SR) 
          qt=quantile(Target_birds[,"SR"], probs=c(0.1, 0.9))
          rownames(Target_birds) <- Target_birds$ID
          Target_birds<-Target_birds[,-1]
          

## Target definition
target_func <- function(SR, qt, log=TRUE){
  SR_i <- SR
  qt_i <- qt
  
  if(log) {
    SR <- log(SR)
    qt <- log(qt)
  }	
  dat <- data.frame(matrix(c(100,10,qt), ncol=2, dimnames=list(NULL, c("Target", "SR"))))
  lm_target <- lm(Target~SR, data=dat)
  Tar <- predict(lm_target,newdata=as.data.frame(SR))
  Tar[SR_i<=qt_i[1]] <- 100
  Tar[SR_i>=qt_i[2]] <- 10
  return(Tar)
}

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

grid.arrange(a,b,ncol=2,top = textGrob("Species target achievements" ,gp=gpar(fontsize=20,font=3)))











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

