#TODO ANNOTE  SCRIPT
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


load(file=file.path(results_dir,"mammals","50km","mammals_PA.RData"))
Target_mammals <- merge (Target_mammals,mammals_PA, by = "row.names" )
rownames(Target_mammals) <- Target_mammals[,1]
Target_mammals <- Target_mammals[,c(2,3,5,6,15)]


Target_mammals[,"TargetExp"] <- target_func(Target_mammals[,"SR"], qt, log=T)
Target_mammals[,"TargetMet_CellsWithPA"] <- 100*(Target_mammals[,"PercentageCellsWithPA"]/Target_mammals[,"TargetExp"])
Target_mammals[,"TargetMet_Percentagecover"] <- 100*(Target_mammals[,"Percentagecover"]/Target_mammals[,"TargetExp"])
Target_mammals <- merge(Target_mammals,data_DR_mammals,by="row.names")
Target_mammals <- na.omit(Target_mammals)

load(file=file.path(results_dir,"birds","50km","birds_PA.RData"))
Target_birds <- merge (Target_birds,birds_PA, by = "row.names" )
rownames(Target_birds) <- Target_birds[,1]
Target_birds <- Target_birds[,c(2,3,5,6,15)]


Target_birds[,"TargetExp"] <- target_func(Target_birds[,"SR"], qt, log=T)
Target_birds[,"TargetMet_CellsWithPA"] <- 100*(Target_birds[,"PercentageCellsWithPA"]/Target_birds[,"TargetExp"])
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
