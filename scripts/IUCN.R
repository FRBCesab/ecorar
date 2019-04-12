#GET THE REDLIST STATUS 

rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NL")

library(rredlist) 
library(ggplot2)
library(gridExtra)
library(ggsignif)

#LOAD DATA ---- 


#----

#REDLIST FUNK RARE ---- 
load(file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))
load(file=file.path(results_dir,"birds","50km","FR_birds.RData"))

load(file=file.path(data_dir,"birds/redlist_birds.RData"))
load(file=file.path(data_dir,"mammals/redlist_mammals.RData"))


load(file=file.path(data_dir,"birds/birdsID.RData"))
load(file=file.path(data_dir,"mammals/mammalsID.RData"))

redlist_birds <- merge(birdsID,redlist_birds,by.x = "Name", by.y = "Scientific")
redlist_birds <- redlist_birds[,c(1,2,4)]
colnames(redlist_birds) <- c("Species","SpecID", "IUCN_status")

redlist_mammals <- merge(mammalsID,redlist_mammals,by.x = "Name", by.y = "Scientific")
#redlist_mammals <- redlist_mammals[,c(4,2,7)]
redlist_mammals <- redlist_mammals[,c(1,2,4)]
colnames(redlist_mammals) <- c("Species","SpecID", "IUCN_status")

plot.iucn <- function(taxa,FR_all,redlist){
  
  #taxa= "birds"
  #FR_all <- FR_birds
  #redlist <- redlist_birds
  
  #taxa= "mammals"
  #FR_all <- FR_mammals
  #redlist <- redlist_mammals
  
  #create the data_iucnFR dataframe
  data_iucnFR <- merge(FR_all$FR,redlist,by.x = "row.names", by.y = "SpecID")
  data_iucnFR <- data_iucnFR[data_iucnFR$IUCN_status %in% c("LC_NT", "NE","THR"),]
  rownames(data_iucnFR) <- data_iucnFR$species
  
  data_iucnFR$DR_class="NA"
  
  QD75 <- FR_all$Q$Q75_D
  QD25 <- FR_all$Q$Q25_D
  QR75 <- FR_all$Q$Q75_R
  QR25 <- FR_all$Q$Q25_R
  
  data_iucnFR$DR_class[(data_iucnFR$Din<QD25) & (data_iucnFR$Rin<QR25)]="D25R25"
  data_iucnFR$DR_class[(data_iucnFR$Din>QD75) & (data_iucnFR$Rin>QR75)]="D75R75"
  data_iucnFR$DR_class[(data_iucnFR$Din<QD25) & (data_iucnFR$Rin>QR75)]="D25R75"
  data_iucnFR$DR_class[(data_iucnFR$Din>QD75) & (data_iucnFR$Rin<QR25)]="D75R25"
  data_iucnFR$DR_class[(((data_iucnFR$Din>QD25) & (data_iucnFR$Din<QD75)) & ((data_iucnFR$Rin>QR25) & (data_iucnFR$Rin<QR75)))]="AVG"
  
  #IUCN & DR status 
  
  
  D25R25 <- as.data.frame(table(data_iucnFR[data_iucnFR$DR_class=="D25R25","IUCN_status"]))
  colnames(D25R25) <- c('IUCN','N')
  D25R25$Freq <- D25R25$N/sum(D25R25$N)
  D25R25 <- cbind(DR="D25R25",D25R25)
 
  D75R75 <- as.data.frame(table(data_iucnFR[data_iucnFR$DR_class=="D75R75","IUCN_status"]))
  colnames(D75R75) <- c('IUCN','N')
  D75R75$Freq <- D75R75$N/sum(D75R75$N)
  D75R75 <- cbind(DR="D75R75",D75R75)
  
  #D25R75 <- as.data.frame(table(data_iucnFR[data_iucnFR$DR_class=="D25R75","IUCN_status"]))
  #colnames(D25R75) <- c('IUCN','N')
  #D25R75$Freq <- D25R75$N/sum(D25R75$N)
  #D25R75 <- cbind(DR="D25R75",D25R75)
  
  #D75R25 <- as.data.frame(table(data_iucnFR[data_iucnFR$DR_class=="D75R25","IUCN_status"]))
  #colnames(D75R25) <- c('IUCN','N')
  #D75R25$Freq <- D75R25$N/sum(D75R25$N)
  #D75R25 <- cbind(DR="D75R25",D75R25)
  
  AVG <- as.data.frame(table(data_iucnFR[data_iucnFR$DR_class=="AVG","IUCN_status"]))
  colnames(AVG) <- c('IUCN','N')
  AVG$Freq <- AVG$N/sum(AVG$N)
  AVG <- cbind(DR="AVG",AVG)
  
  #DR_IUCN <- rbind(D25R25,AVG,D75R25,D25R75,D75R75)
  DR_IUCN <- rbind(AVG,D25R25,D75R75)
  if(taxa=="mammals")  DR_IUCN <- rbind(DR_IUCN,data.frame(DR="D25R25",IUCN="NE",N="NA",Freq=0))
  #Plot#1 
  size=0.2
  width=0.2
  
  #a <- ggplot(data_iucnFR, aes(x=IUCN_status, y=Rin, fill=IUCN_status)) + geom_boxplot() + guides(fill=FALSE) +
  # geom_jitter(width = width,size=size,color="darkgrey") + scale_y_continuous(limits = c(0.8, 1.15)) + 
  #  ggsignif::geom_signif(comparisons = list(c("LC_NT", "NE"),c("NE", "THR"),c("LC_NT", "THR")),y_position = c(1.05,1.1,1.15),map_signif_level=TRUE,tip_length=0.01)
  
  # b <- ggplot(data_iucnFR, aes(x=IUCN_status, y=Din, fill=IUCN_status)) + geom_boxplot() + guides(fill=FALSE) +
  #  geom_jitter(width = width,size=size,color="darkgrey") + scale_y_continuous(limits = c(0.15, 1.15)) + 
  #  ggsignif::geom_signif(comparisons = list(c("LC_NT", "NE"),c("NE", "THR"),c("LC_NT", "THR")),y_position = c(1.05,1.1,1.15),map_signif_level=TRUE,tip_length=0.01)

  # pdf(file.path(results_dir,paste0(taxa,"/","50km","/figs/IUCN#1.pdf")),width=12,height=4)
  #grid.arrange(a,b,ncol=2)
  # dev.off()
  
  #Plot#2
  DR_IUCN$IUCN<- factor(DR_IUCN$IUCN,levels = c('LC_NT','NE','THR'),ordered = TRUE)
  #pdf(file.path(results_dir,paste0(taxa,"/","50km","/figs/IUCN#2.pdf")),width=6,height=4)
  ggplot(data=DR_IUCN, aes(x=IUCN, y=Freq, fill=DR)) + scale_fill_manual(values=c("#00AFBB", "#E7B800","orangered"))+#    scale_fill_manual(values=c("#5E4FA2",  "#66C2A5",  "#FDAE61", "#D53E4F", "#9E0142")) + 
  geom_bar(stat="identity", position=position_dodge()) + ggtitle("IUCN Status") + theme_bw() + scale_y_continuous(name="Frequence", limits=c(0, 1))
  #dev.off() 
}




plot.iucn(taxa="birds",FR_all=FR_birds,redlist=redlist_birds)
plot.iucn(taxa="mammals",FR_all=FR_mammals,redlist=redlist_mammals)

