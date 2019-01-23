#PCO_VISU

library(ggplot2)
library(ape)
library(dplyr)
library(gridExtra)  


rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NM")


#LOAD DATA ---- 

##Birds 

load(file=file.path(results_dir,"birds/50km/disTraits_birds.RData")) #TODO why the data are in 50km ? 
load(file=file.path(results_dir,"birds//50km/birdstrait.RData")) #TODO why the data are in 50km ? 
load(file=file.path(results_dir,"birds/50km/FR_birds.RData"))

##Mammals 

load(file=file.path(results_dir,"mammals/50km/disTraits_mammals.RData")) #TODO why the data are in 50km ? 
load(file=file.path(results_dir,"mammals/50km/mammalstrait.RData")) #TODO why the data are in 50km ? 
load(file=file.path(results_dir,"mammals/50km/FR_mammals.RData"))

#----

#COMPUTE PCOA ---- 

  ##Birds
    pco_birds<-pcoa(disTraits_birds)
    save(pco_birds, file=file.path(results_dir,"birds/pco_birds.RData"))
  
  ##Mammals 
    pco_mammals<-pcoa(disTraits_mammals)
    save(pco_mammals, file=file.path(results_dir,"mammals/pco_mammals.RData"))

#----
    
    
#LOAD PCOA ---- 
    
  load(file=file.path(results_dir,"birds/pco_birds.RData"))
  load(file=file.path(results_dir,"mamals/pco_mamals.RData"))
  
#----
  
#INTERPRETING THE PCOA ----  
  
  pco_data <- pco_birds
  traits <- birdstrait
  FR_data <- FR_birds
  
  ###histo of first eigenvalues 
  
    df <- data.frame(Eigenvalues = paste0('E',seq(1:8)),
                     Value = pco_data$values$Relative_eig[1:8])
    
    ggplot(data=df, aes(x=Eigenvalues, y=Value)) + geom_bar(stat="identity",fill="steelblue")+
      geom_text(aes(label=round(Value,digits=2)), vjust=1.6, color="white", size=3.5)+
      theme_minimal()
  
  ###plot pcoa
  
    jitval=500
    
    df <- data.frame(x = jitter(pco_data$vectors[,"Axis.1"],jitval),
                     y = jitter(pco_data$vectors[,"Axis.2"],jitval))
    pco12 <- ggplot(df, aes(x, y)) + geom_point(colour = "black")+
      labs(x = "PC1",y = "PC2")+ theme_minimal() +
      geom_hline(yintercept=0,linetype="dashed", color = "red", size=0.5) + 
      geom_vline(xintercept=0,linetype="dashed", color = "red", size=0.5) 
    
    df <- data.frame(x = jitter(pco_data$vectors[,"Axis.1"],jitval),
                     y = jitter(pco_data$vectors[,"Axis.3"],jitval))
    pco13 <- ggplot(df, aes(x, y)) + geom_point(colour = "black")+
      labs(x = "PC1",y = "PC3")+ theme_minimal()+
      geom_hline(yintercept=0,linetype="dashed", color = "red", size=0.5) + 
      geom_vline(xintercept=0,linetype="dashed", color = "red", size=0.5) 
    
    df <- data.frame(x = jitter(pco_data$vectors[,"Axis.2"],jitval),
                     y = jitter(pco_data$vectors[,"Axis.3"],jitval))
    pco23 <- ggplot(df, aes(x, y)) + geom_point(colour = "black")+
      labs(x = "PC2",y = "PC3")+ theme_minimal()+
      geom_hline(yintercept=0,linetype="dashed", color = "red", size=0.5) + 
      geom_vline(xintercept=0,linetype="dashed", color = "red", size=0.5) 
    
    grid.arrange(pco12,pco13,pco23, nrow=2,ncol=2)
  
  ###Corelations with the traits 
  
    trait_table = traits %>% mutate_all(as.numeric)
    traits_pcoa <- merge(traits,pco_data$vectors,by="row.names",all.x=TRUE)
    
    cor_test <- do.call(rbind,lapply(colnames(traits),fun <- function(vname){
      #vname <- colnames(mamalstrait)[1]
      
      if(is.numeric(traits_pcoa[1,vname])){
        
        trait.st = apply(trait_table, 2, scale, center=TRUE, scale=TRUE)
        ax1 <- cor(scale(traits_pcoa[,vname],center=TRUE, scale=TRUE),traits_pcoa$Axis.1, method="spearman")
        ax2 <- cor(scale(traits_pcoa[,vname],center=TRUE, scale=TRUE),traits_pcoa$Axis.2, method="spearman")
        ax3 <- cor(scale(traits_pcoa[,vname],center=TRUE, scale=TRUE),traits_pcoa$Axis.3, method="spearman")
        test <- "spearman"
      } else {
        
        nfree <- length(unique(traits_pcoa[,vname]))-1
        
        ktest <- kruskal.test(traits_pcoa$Axis.1 ~ traits_pcoa[,vname])
        Fval <- as.numeric(ktest$statistic)/nfree
        ax1 <- (Fval*nfree)/(Fval*nfree+dim(traits)[1]-nfree)
        
        ktest <- kruskal.test(traits_pcoa$Axis.2 ~ traits_pcoa[,vname])
        Fval <- as.numeric(ktest$statistic)/nfree
        ax2 <- (Fval*nfree)/(Fval*nfree+dim(traits)[1]-nfree)
        
        ktest <- kruskal.test(traits_pcoa$Axis.3 ~ traits_pcoa[,vname])
        Fval <- as.numeric(ktest$statistic)/nfree
        ax3 <- (Fval*nfree)/(Fval*nfree+dim(traits)[1]-nfree)
        test <- "Kruskal"
        
      }
      
      cbind.data.frame(var=vname,test=test ,ax1=ax1 ,ax2=ax2,ax3=ax3)
      
    }))
  
  
  ###Plot all TO DO A reprendre pour que cela ne soit pas trop degeu ... 
  
    ttheme_minimal(base_size = 12, base_colour = "black", base_family = "",
                   parse = FALSE, padding = unit(c(4, 4), "mm"))
    tbl <- tableGrob(cor_test, rows=NULL, theme = ttheme_minimal(base_size = 08))
    
    grid.arrange(hist_eigen, tbl,pco12,pco13,
                 nrow=2,ncol=2,
                 as.table=TRUE)
  
  ##Correlation between traits and pcoa vs. distinctivness & restrictedness TO DO cor_test_FR ne marche pas ... 
  
    var="Din"
    FR_pcoa <- merge(FR_data$FR,pco_data$vectors,by="row.names",all.x=TRUE)
  
  a <- ggplot(FR_pcoa, aes(x=Axis.1, y=get(var))) + geom_point() + stat_smooth(method = "lm", size = 1) + labs(x = "PC1",y = var)+ theme_minimal()
  b <- ggplot(FR_pcoa, aes(x=Axis.2, y=get(var))) + geom_point() + stat_smooth(method = "lm", size = 1) + labs(x = "PC2",y = var)+ theme_minimal()
  c <- ggplot(FR_pcoa, aes(x=Axis.3, y=get(var))) + geom_point() + stat_smooth(method = "lm", size = 1) + labs(x = "PC3",y = var)+ theme_minimal()
  
  traits_pcoa <- merge(traits,pco_data$vectors,by="row.names",all.x=TRUE)
  
  cor_test_FR <- do.call(rbind,lapply(colnames(traits),fun <- function(vname,var="Din"){
    
    
    ##var="Din"
    ##vname <- colnames(mamalstrait)[1]
    
    traits_FR <- merge(traits,FR_data$FR,by="row.names",all.x=TRUE)
    
    if(is.numeric(traits_FR[1,vname])){
      
      Din <- cor(scale(traits_FR[,vname],center=TRUE, scale=TRUE),traits_FR$Din, method="spearman")
      Rin <- cor(scale(traits_FR[,vname],center=TRUE, scale=TRUE),traits_FR$Rin, method="spearman")
      test <- "spearman"
    } else {
      
      nfree <- length(unique(traits_FR[,vname]))-1
      
      ktest <- kruskal.test(traits_FR$Din ~ traits_FR[,vname])
      Fval <- as.numeric(ktest$statistic)/nfree
      Din <- (Fval*nfree)/(Fval*nfree+dim(mamalstrait)[1]-nfree)
      
      ktest <- kruskal.test(traits_FR$Rin ~ traits_FR[,vname])
      Fval <- as.numeric(ktest$statistic)/nfree
      Rin <- (Fval*nfree)/(Fval*nfree+dim(mamalstrait)[1]-nfree)
      
      test <- "Kruskal"
      
    }
    
    cbind.data.frame(var=vname,test=test ,Din=Din ,Rin=Rin)
    
  }))
  
  ttheme_minimal(base_size = 12, base_colour = "red", base_family = "",
                 parse = FALSE, padding = unit(c(2, 2), "mm"))
  tbl <- tableGrob(cor_test_FR, rows=NULL, theme = ttheme_minimal(base_size = 04))
  
  grid.table(cor_test_FR)
  
  grid.arrange(a, b,c,ncol=1)