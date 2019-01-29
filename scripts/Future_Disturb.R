#Future
library(raster)
library(ggplot2)
library(gridExtra)
library(grid)
library(foreign)
library(ggsignif)
library(RColorBrewer)
library(viridis)
#

#LOST OF CELLS CLIMATE CHANGE --- 

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
    geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(-c, ymax)) + geom_hline(yintercept=0,col="red",linetype="dashed") + 
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

#----



#LAND USES --- 
#----
#LOAD & FORMAT DATA --- 
deltaShortVeg.tif
deltaCanopy.tif
deltaBare.tif

deltaShortVeg <- readGDAL(file=file.path(data_dir,"Environmental",reso,"deltaShortVeg.tif"))
raster("a.tif")

rgb <- brick(“pathto/rgb.tif”)

#----

# PROTECTED----
#LOAD & FORMAT DATA --- 
load(file=file.path(data_dir,"Environmental",reso,"PA&Country.RData"))
load(file=file.path(results_dir,"birds","50km","funk_birds.RData"))

# Pour chaque espèce: Il faut le nombre cells protégées, nombre cell tot, le pourcentage de couverture de no take, le nombre de conflit moyen, 
species<-unique(unlist(occ_birds_list))
birds_PA<- do.call(rbind,mclapply(species, function(id) {
  cells <- sapply(1:length(occ_birds_list), function(i) any(occ_birds_list[[i]] == id))
  names(cells)<-names(occ_birds_list)
  cells<-PA_Country[PA_Country$ID %in% names(cells[cells=="TRUE"]),]
  n_cell<-dim(cells)[1]
  n_pro<-sum(na.omit(cells$PA_DEF==1))
  PA_mean<-apply(cells[,c(30,32,34)],2,mean,na.rm=T)
  info<-t(data.frame(c(n_cell,n_pro,PA_mean)))
  colnames(info)<-c("N_cells","N_cells_protected","Percentagecover","meanConflict","meanHDI")
  rownames(info)<-id
  return(info)
},mc.cores=3))

birds_PA<-merge(birds_PA,FR_birds$FR,by="row.names")
rownames(birds_PA) <- birds_PA[,1]
birds_PA <- birds_PA[,-1]
birds_PA$PercentageCellsWithPA<- (birds_PA$N_cells_protected/birds_PA$N_cells)*100
save(birds_PA,file=file.path(results_dir,"birds",reso,"birds_PA.RData")) 


species<-unique(unlist(occ_mammals_list))
mammals_PA<- do.call(rbind,mclapply(species, function(id) {
  cells <- sapply(1:length(occ_mammals_list), function(i) any(occ_mammals_list[[i]] == id))
  names(cells)<-names(occ_mammals_list)
  cells<-PA_Country[PA_Country$ID %in% names(cells[cells=="TRUE"]),]
  n_cell<-dim(cells)[1]
  n_pro<-sum(na.omit(cells$PA_DEF==1))
  PA_mean<-apply(cells[,c(30,32,34)],2,mean,na.rm=T)
  info<-t(data.frame(c(n_cell,n_pro,PA_mean)))
  colnames(info)<-c("N_cells","N_cells_protected","Percentagecover","meanConflict","meanHDI")
  rownames(info)<-id
  return(info)
},mc.cores=3))

mammals_PA<-merge(mammals_PA,FR_mammals$FR,by="row.names")
rownames(mammals_PA) <- mammals_PA[,1]
mammals_PA <- mammals_PA[,-1]
mammals_PA$PercentageCellsWithPA<- (mammals_PA$N_cells_protected/mammals_PA$N_cells)*100

save(mammals_PA,file=file.path(results_dir,"mammals",reso,"mammals_PA.RData")) 

plot_PA <- function(taxa,FR_all,data_PA){  
   #taxa="mammals"
   #FR_all=FR_mammals
   #data_PA=mammals_PA

  data_PA$DR_class="NA"
  
  QD75 <- FR_all$Q$Q75_D
  QD25 <- FR_all$Q$Q25_D
  QR75 <- FR_all$Q$Q75_R
  QR25 <- FR_all$Q$Q25_R
  
  data_PA$DR_class[(data_PA$Din<QD25) & (data_PA$Rin<QR25)]="D25R25"
  data_PA$DR_class[(data_PA$Din>QD75) & (data_PA$Rin>QR75)]="D75R75"
  data_PA$DR_class[(data_PA$Din<QD25) & (data_PA$Rin>QR75)]="D25R75"
  data_PA$DR_class[(data_PA$Din>QD75) & (data_PA$Rin<QR25)]="D75R25"
  data_PA$DR_class[(((data_PA$Din>QD25) & (data_PA$Din<QD75)) & ((data_PA$Rin>QR25) & (data_PA$Rin<QR75)))]="AVG"

  data_PA$InvRin=1-data_PA$Rin
  
  data_plot <- data_PA[data_PA$DR_class!='NA',]
  
  #col_br <- brewer.pal(n = 5, name = "Spectral")
  
  col_br <- viridis(n = 6, option = "A") # for colorblind people
  
  ymax=100
  a <- ggplot(data_plot, aes(x=DR_class, y=Percentagecover, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
    geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(data_plot$Percentagecover,na.rm=T),col="red",linetype="dashed") + 
    ggsignif::geom_signif(comparisons = list(c("D25R25", "D75R75"),c("D25R25", "D25R75"),c("D25R25", "D75R25")),y_position = c(ymax-25,ymax-10,ymax),map_signif_level=TRUE,tip_length=0.01)+
    labs(x = "DR class",y="Percentage")
  
  data_plot_sub <- data_plot[((data_plot$DR_class=='D25R25') | (data_plot$DR_class=='D75R75') | (data_plot$DR_class=='AVG')),]
  b <- ggplot(data_plot_sub, aes(Percentagecover,fill=DR_class,color=DR_class)) + geom_density(adjust = 1.5,alpha = 0.1) + xlim(0, ymax)+ 
    scale_fill_manual(values=c(col_br[1],col_br[2],col_br[5]))+ scale_color_manual(values=c(col_br[1],col_br[2],col_br[5]))+
    theme(legend.position = c(0.9, 0.8)) + geom_vline(xintercept=mean(data_plot_sub$Percentagecover,na.rm=T),col="red",linetype="dashed")+
    labs(x = "Percentage")
 
  pdf(file.path(results_dir,paste0(taxa,"/50km","/figs/Average cover of protected area per cells.pdf")),width=12,height=8) 
  grid.arrange(a,b,ncol=2,top = textGrob("Average cover of protected area per cells" ,gp=gpar(fontsize=20,font=3)))
  dev.off()

  
  ymax=1.1
  c <- ggplot(data_plot, aes(x=DR_class, y=meanHDI, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
    geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(data_plot$meanHDI,na.rm=T),col="red",linetype="dashed") + 
    ggsignif::geom_signif(comparisons = list(c("D25R25", "D75R75"),c("D25R25", "D25R75"),c("D25R25", "D75R25")),y_position = c(ymax-0.1,ymax-0.05,ymax),map_signif_level=TRUE,tip_length=0.01)+
    labs(x = "DR class",y="HDI")
  
  data_plot_sub <- data_plot[((data_plot$DR_class=='D25R25') | (data_plot$DR_class=='D75R75') | (data_plot$DR_class=='AVG')),]
  d <- ggplot(data_plot_sub, aes(meanHDI,fill=DR_class,color=DR_class)) + geom_density(adjust = 1.5,alpha = 0.1) + xlim(0, ymax)+ 
    scale_fill_manual(values=c(col_br[1],col_br[2],col_br[5]))+ scale_color_manual(values=c(col_br[1],col_br[2],col_br[5]))+
    theme(legend.position = c(0.9, 0.8)) + geom_vline(xintercept=mean(data_plot_sub$meanHDI,na.rm=T),col="red",linetype="dashed")+
    labs(x = "HDI")
  
  pdf(file.path(results_dir,paste0(taxa,"/50km","/figs/HDI.pdf")),width=12,height=8) 
  grid.arrange(c,d,ncol=2,top = textGrob("Human Development Index" ,gp=gpar(fontsize=20,font=3)))
  dev.off()
  
  
  ymax=260
  e <- ggplot(data_plot, aes(x=DR_class, y=meanConflict, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
    geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(data_plot$meanConflict,na.rm=T),col="red",linetype="dashed") + 
    ggsignif::geom_signif(comparisons = list(c("D25R25", "D75R75"),c("D25R25", "D25R75"),c("D25R25", "D75R25")),y_position = c(ymax-25,ymax-12,ymax),map_signif_level=TRUE,tip_length=0.01)+
    labs(x = "DR class",y="Number of conflicts")
  
  data_plot_sub <- data_plot[((data_plot$DR_class=='D25R25') | (data_plot$DR_class=='D75R75') | (data_plot$DR_class=='AVG')),]
  f <- ggplot(data_plot_sub, aes(meanConflict,fill=DR_class,color=DR_class)) + geom_density(adjust = 1.5,alpha = 0.1) + xlim(0, ymax)+ 
    scale_fill_manual(values=c(col_br[1],col_br[2],col_br[5]))+ scale_color_manual(values=c(col_br[1],col_br[2],col_br[5]))+
    theme(legend.position = c(0.9, 0.8)) + geom_vline(xintercept=mean(data_plot_sub$meanConflict,na.rm=T),col="red",linetype="dashed")+
    labs(x = "Number of conflicts")
  
  pdf(file.path(results_dir,paste0(taxa,"/50km","/figs/Conflict.pdf")),width=12,height=8) 
  grid.arrange(e,f,ncol=2,top = textGrob("Average number of conflicts" ,gp=gpar(fontsize=20,font=3)))
  dev.off()
  
  ymax=150
  g <- ggplot(data_plot, aes(x=DR_class, y=PercentageCellsWithPA, fill=DR_class)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_manual(values=col_br)+
    geom_jitter(width = 0.1,size=0.5,color="darkgrey") + scale_y_continuous(limits = c(0, ymax)) + geom_hline(yintercept=mean(data_plot$PercentageCellsWithPA,na.rm=T),col="red",linetype="dashed") + 
    ggsignif::geom_signif(comparisons = list(c("D25R25", "D75R75"),c("D25R25", "D25R75"),c("D25R25", "D75R25")),y_position = c(ymax-25,ymax-12,ymax),map_signif_level=TRUE,tip_length=0.01)+
    labs(x = "DR class",y="Number of conflicts")
  
  data_plot_sub <- data_plot[((data_plot$DR_class=='D25R25') | (data_plot$DR_class=='D75R75') | (data_plot$DR_class=='AVG')),]
  h <- ggplot(data_plot_sub, aes(PercentageCellsWithPA,fill=DR_class,color=DR_class)) + geom_density(adjust = 1.5,alpha = 0.1) + xlim(0, ymax)+ 
    scale_fill_manual(values=c(col_br[1],col_br[2],col_br[5]))+ scale_color_manual(values=c(col_br[1],col_br[2],col_br[5]))+
    theme(legend.position = c(0.9, 0.8)) + geom_vline(xintercept=mean(data_plot_sub$PercentageCellsWithPA,na.rm=T),col="red",linetype="dashed")+
    labs(x = "Number of conflicts")
  
  pdf(file.path(results_dir,paste0(taxa,"/50km","/figs/PercentageCellsWithPA.pdf")),width=12,height=8) 
  grid.arrange(g,h,ncol=2,top = textGrob("Percentage of Cells With PA" ,gp=gpar(fontsize=20,font=3)))
  dev.off()
  
  }
  

plot_PA(taxa="birds",FR_birds,birds_PA)
plot_PA(taxa="mammals",FR_mammals,mammals_PA)

