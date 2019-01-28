#Future

library(ggplot2)
library(gridExtra)
library(grid)
library(foreign)
library(ggsignif)
library(RColorBrewer)
library(viridis)



#LOAD & FORMAT DATA ---- 

##Mammals
taxocor <- read.csv2(file.path(data_dir,"mammals","TaxonomicCorrespondancesMammals.csv"),header=TRUE)

iucn_code <- data.frame(Scientific=taxocor$DistriName, IUCN_code=taxocor$DistriCode)
iucn_code <- merge(mammalsID,iucn_code,by.x="Name", by.y="Scientific")
rownames(iucn_code) <-iucn_code$IUCN_code 

mammals_future <- read.table(file.path(data_dir,"mammals","ENSEMBLE_CA_Mammals_LossStableGain_CurrFut.txt"),header=TRUE)
scenar <- unique(mammals_future$SCE)

mammals_future_scenar_all <- lapply(scenar, function(id) {
  #id=scenar[1]
  mammals_future_scenar <- subset(mammals_future,mammals_future$SCE==id)
  mammals_future_scenar_mean <- aggregate(mammals_future_scenar[, 6:10], list(mammals_future_scenar$SP), mean)
  rownames(mammals_future_scenar_mean) <- mammals_future_scenar_mean$Group.1
  mammals_future_scenar_mean <- mammals_future_scenar_mean[,-1]
  mammals_future_scenar_all <- merge(iucn_code,mammals_future_scenar_mean,by="row.names")
  rownames(mammals_future_scenar_all) <- mammals_future_scenar_all$ID
  mammals_future_scenar_all[,-1]
})

names(mammals_future_scenar_all) <- scenar

load(file=file.path(results_dir,"mammals","50km","FR_mammals.RData"))


##Birds

load(file=file.path(results_dir,"birds",reso,"FR_birds.RData"))

birds_future <- read.table(file.path(data_dir,"birds","ENSEMBLE_CA_Birds_LossStableGain_CurrFut.txt"),header=TRUE)
birds_future$SP <- paste0("sp",birds_future$SP)

scenar <- unique(birds_future$SCE)

birds_future_scenar_all <- lapply(scenar, function(id) {
  #id=scenar[1]
  birds_future_scenar <- subset(birds_future,birds_future$SCE==id)
  birds_future_scenar_mean <- aggregate(birds_future_scenar[, 6:10], list(birds_future_scenar$SP), mean)
  rownames(birds_future_scenar_mean) <- birds_future_scenar_mean$Group.1
  birds_future_scenar_mean <- birds_future_scenar_mean[,-1]
  birds_future_scenar_all <- merge(birdsID,birds_future_scenar_mean,by.x="ID", by.y="row.names")
  rownames(birds_future_scenar_all) <- birds_future_scenar_all$ID
  birds_future_scenar_all[,-1]
})

names(birds_future_scenar_all) <- scenar
colnames(birds_future_scenar_all)[1]  <- 
load(file=file.path(results_dir,"birds/FR_birds_all.RData"))


#----

#FUTUR ----

plot_futur <- function(taxa,FR_all,id_scenar,futur_all,ymax)
{
  
  # taxa="mammals"
  # FR_all=FR_mammals
  # id_scenar=scenar[11]
  # futur_all=mammals_future_scenar_all
  # ymax <- 300

  # taxa="birds"
  # FR_all=FR_birds
  # id_scenar=scenar[11]
  # futur_all=birds_future_scenar_all
  # ymax <- 300
  
  future <- futur_all[[id_scenar]]
  #rownames(future) <- future$ID
  
  data_future <- merge(FR_all$FR,future,by="row.names")
  rownames(data_future) <- data_future$Row.names
  
  
  data_future$delta=100*(data_future$FUTUR-data_future$CURRENT)/data_future$CURRENT
  
  data_future$DR_class="NA"
  
  QD75 <- FR_all$Q$Q75_D
  QD25 <- FR_all$Q$Q25_D
  QR75 <- FR_all$Q$Q75_R
  QR25 <- FR_all$Q$Q25_R
  
  data_future$DR_class[(data_future$Din<QD25) & (data_future$Rin<QR25)]="D25R25"
  data_future$DR_class[(data_future$Din>QD75) & (data_future$Rin>QR75)]="D75R75"
  data_future$DR_class[(data_future$Din<QD25) & (data_future$Rin>QR75)]="D25R75"
  data_future$DR_class[(data_future$Din>QD75) & (data_future$Rin<QR25)]="D75R25"
  data_future$DR_class[(((data_future$Din>QD25) & (data_future$Din<QD75)) & ((data_future$Rin>QR25) & (data_future$Rin<QR75)))]="AVG"
  
  #sum(data_future$DR_class=="D25R25")
  #sum(data_future$DR_class=="D75R75")
  #sum(data_future$DR_class=="D25R75")
  #sum(data_future$DR_class=="D75R25")
  #sum(data_future$DR_class=="AVG")
  
  data_future$InvRin=1-data_future$Rin
  
  #Plot#1
  a <- ggplot(data_future, aes(x=CURRENT, y=delta)) + geom_point() + stat_smooth(method = "lm", size = 1) + geom_hline(yintercept=0,col="red",linetype="dashed") +scale_x_continuous(trans='log10')+ scale_y_continuous(limits = c(-100, ymax))
  
  b <- ggplot(data_future, aes(x=Din, y=delta)) + geom_point() + stat_smooth(method = "lm", size = 1) + geom_hline(yintercept=0,col="red",linetype="dashed")+ scale_y_continuous(limits = c(-100, ymax))
  
  
  #Plot#2
  
  data_plot <- data_future[data_future$DR_class!='NA',]
  
  #col_br <- brewer.pal(n = 5, name = "Spectral")
  
  col_br <- viridis(n = 6, option = "A") # for colorblind people
  
  c <- ggplot(data_plot, aes(x=DR_class, y=delta, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
    geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(-100, ymax)) + geom_hline(yintercept=0,col="red",linetype="dashed") + 
    ggsignif::geom_signif(comparisons = list(c("D25R25", "D75R75"),c("D25R25", "D25R75"),c("D25R25", "D75R25")),y_position = c(ymax-100,ymax-50,ymax),map_signif_level=TRUE,tip_length=0.01)
  
  data_plot_sub <- data_plot[((data_plot$DR_class=='D25R25') | (data_plot$DR_class=='D75R75') | (data_plot$DR_class=='AVG')),]
  d <- ggplot(data_plot_sub, aes(delta,fill=DR_class,color=DR_class)) + geom_density(adjust = 1.5,alpha = 0.1) + xlim(-100, ymax)+ 
    scale_fill_manual(values=c(col_br[1],col_br[2],col_br[5]))+ scale_color_manual(values=c(col_br[1],col_br[2],col_br[5]))+
    theme(legend.position = c(0.9, 0.8)) + geom_vline(xintercept=0,col="red",linetype="dashed") 
  
  grid.arrange(a,b,c,d,ncol=2,top = textGrob(id_scenar,gp=gpar(fontsize=20,font=3)))
  
}

plot_futur(taxa="birds",FR_all=FR_birds,id_scenar=scenar[4],futur_all=birds_future_scenar_all,ymax= 300)

plot_futur(taxa="mammals",FR_all=FR_mammals,id_scenar=scenar[4],futur_all=mammals_future_scenar_all,ymax= 300)


pdf(file.path(results_dir,paste0(taxa,"/50km","/figs/FUTUR.pdf")),width=12,height=8)
for (i in 1:length(scenar)){
  plot_futur(taxa="birds",FR_all=FR_birds,id_scenar=scenar[i],futur_all=birds_future_scenar_all,ymax= 300)
}
dev.off()  

taxa="mamals"
pdf(file.path(results_dir,paste0(taxa,"/figs/FUTUR.pdf")),width=12,height=8)
for (i in 1:length(scenar)){
  plot_futur(taxa=taxa,FR_all=FR_mamals_all,id_scenar=scenar[i],futur_all=mammals_future_scenar_all,ymax= 300)
}
dev.off()  







