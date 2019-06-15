#PCO_VISU

rm(list=ls(all=TRUE)) 


library(ggplot2)
library(ape)
library(dplyr)
library(gridExtra)  
library(plyr)

source("./scripts/Functions.R")
who.remote(remote=TRUE,who="NM")


#LOAD DATA ---- 

##Birds 

load(file=file.path(results_dir,"birds/disTraits_birds.RData")) 
load(file=file.path(results_dir,"birds/birdstrait.RData"))
load(file=file.path(results_dir,"birds/50km/FR_birds.RData"))

##Mammals 

load(file=file.path(results_dir,"mammals/disTraits_mammals.RData")) 
load(file=file.path(results_dir,"mammals/mammalstrait.RData"))
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
  load(file=file.path(results_dir,"mammals/pco_mammals.RData"))
  
#----
  
#INTERPRETING THE PCOA ----  
  
  ###Chosing between birds and mammals
    
    # pco_data <- pco_birds
    # traits <- birdstrait
    # FR_data <- FR_birds
    # taxa <- "birds"
    
    pco_data <- pco_mammals
    traits <- mammalstrait
    FR_data <- FR_mammals
    taxa <- "mammals"
  
  ###histo of first eigenvalues 
  
    df <- data.frame(Eigenvalues = paste0('E',seq(1:8)),
                     Value = pco_data$values$Relative_eig[1:8])
    
    hist_eigen <- ggplot(data=df, aes(x=Eigenvalues, y=Value)) + geom_bar(stat="identity",fill="steelblue")+
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
    
    grid.arrange(hist_eigen,pco12,pco13,pco23, nrow=2,ncol=2)
  
  ###Corelations between traits and pcoa axes 
  
    #TODO verifier type de donnée numeric et al... 
    trait_table = traits %>% mutate_all(as.numeric)
    traits_pcoa <- merge(traits,pco_data$vectors,by="row.names",all.x=TRUE)
    
    cor_test <- do.call(rbind,lapply(colnames(traits),fun <- function(vname){
      
      #vname <- colnames(traits)[1]
      
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
   
    dev.off()
    grid.table(cor_test)
    
  ##Correlation between pcoa axes and Di or Ri 
    
    cor_pcoa <- function(data,pco,var){
      # data=FR_data
      # pco=pco_data
      # var="Din"
      
      df <- data.frame(x1 = pco$vectors[,1],
                       x2 = pco$vectors[,2],
                       x3 = pco$vectors[,3],
                       x4 = pco$vectors[,4],
                       z = data$FR[rownames(pco$vectors),var])
      
      a <- ggplot(df, aes(x=x1, y=z)) + geom_point() + stat_smooth(method = "lm", size = 1) + labs(x = "PC1",y = var)+ theme_minimal()
      b <- ggplot(df, aes(x=x2, y=z)) + geom_point() + stat_smooth(method = "lm", size = 1) + labs(x = "PC2",y = var)+ theme_minimal()
      c <- ggplot(df, aes(x=x3, y=z)) + geom_point() + stat_smooth(method = "lm", size = 1) + labs(x = "PC3",y = var)+ theme_minimal()
      d <- ggplot(df, aes(x=x4, y=z)) + geom_point() + stat_smooth(method = "lm", size = 1) + labs(x = "PC4",y = var)+ theme_minimal()
      
      
      grid.arrange(a,b,c,d,ncol=2)
    }
    
    cor_pcoa(data=FR_data,pco=pco_data,var="Din")
    
    cor_pcoa(data=FR_data,pco=pco_data,var="Rin")
  
  
  ###Exemple of species #TODO A reprendre .... 
    
    sub <- subset(traits_pcoa,(traits_pcoa$Axis.2 < -0.15 & traits_pcoa$Axis.3 < -0.15))
    
    #sub <- subset(traits_pcoa,((traits_pcoa$Axis.2 < -0.26) & ((traits_pcoa$Axis.3 > -0.02) & (traits_pcoa$Axis.3 < 0.02))))
    
    #sub <- subset(traits_pcoa,((traits_pcoa$Axis.3 < -0.25) & ((traits_pcoa$Axis.2 > -0.02) & (traits_pcoa$Axis.2 < 0.02))))
    
    idtotaxo(id=sample(sub$Row.names,1), who="mammals") # TODO la fonction idtotaxo ne marche plus (attention aussi au 50km pour les birds)
    
    ###subspecies 
    
    subspecies<-function(data,pco,axis.x,axis.y,var1,var2,Q1,Q2,DR){
      
      #data=FR_data
      #pco=pco_data
      #axis.x=1
      #axis.y=2
      #var1="Din"
      #var2="Rin"
      #Q1="Q75_D"
      #Q2="Q75_R"
      #DR="DR75"
      
      df <- data.frame(x = pco$vectors[,axis.x],
                       y = pco$vectors[,axis.y],
                       z1 = data$FR[rownames(pco$vectors),var1],
                       z2 = data$FR[rownames(pco$vectors),var2])
      
      if (DR=="DR75") df <- data.frame(df,w=df$z1>as.numeric(data$Q[Q1][[1]]) & df$z2>as.numeric(data$Q[Q2][[1]]))
      if (DR=="DR25") df <- data.frame(df,w=df$z1<as.numeric(data$Q[Q1][[1]]) & df$z2<as.numeric(data$Q[Q2][[1]]))
      
      names(df)[5]<-"w"
      
      df[df$w==TRUE,] 
      
    }
    
    df_sub <- subspecies(data=FR_data, pco=pco_data, axis.x=2, axis.y=3,
                         var1="Din", var2="Rin", Q1="Q25_D", Q2="Q25_R", DR="DR25")
    
    sub <- subset(df_sub,(df_sub$x < 0 & df_sub$y < -0.15))
    
    sub <- subset(df_sub,(((df_sub$x < -0.1) & (df_sub$x > -0.2)) & ((df_sub$y < -0.1) & (df_sub$y > -0.14))))
    
    idtotaxo(id=sample(rownames(sub),1), who="mammals") # TODO la fonction idtotaxo ne marche plus 
    
    head(df_sub)
#----
      
#VISU----
    
  ##Chosing between birds and mammals
    
     # pco_data <- pco_birds
     # traits <- birdstrait
     # FR_data <- FR_birds
     # taxa <- "birds"
    
    pco_data <- pco_mammals
    traits <- mammalstrait
    FR_data <- FR_mammals
    taxa <- "mammals"
    

  ##Plot pcoa for Din & Rin; the function can only be called for var with quantiles  
    
    pcoa.funk.onevar<-function(data,pco,var,plotpdf,resultdir,axis.x,axis.y,jitval,pts,hull,Q1){
      
      # data=FR_data
      # pco=pco_data
      # var="Din"
      # plotpdf=FALSE
      # resultdir=taxa
      # axis.x=1
      # axis.y=2
      # jitval=500
      # pts=TRUE
      # hull=TRUE
      # Q1="Q90_D"
      
      df <- data.frame(x = jitter(pco$vectors[,axis.x],jitval),
                       y = jitter(pco$vectors[,axis.y],jitval),
                       z = data$FR[rownames(pco$vectors),var])
      
      quant <-as.numeric(data$Q[Q1][[1]]) 
      df <- data.frame(df,w=df$z>quant)
      names(df)[4]<-"w"
      df2 <- df[df$z>quant,]

      df2 <- df2[complete.cases(df2), ] # needed because there is one NA in the birds dataframe
      
      
      find_hull <- function(df2) df2[chull(df2$x, df2$y), ]
      hulls <- ddply(df2, "w", find_hull)
      
      p <- ggplot(df, aes(x, y)) +
        geom_point(aes(colour = df$z))+
        scale_colour_gradientn(colours=c("blue","green", "red"),name=var) +
        labs(x = paste0("PC",axis.x),y = paste0("PC",axis.y))+ theme_minimal() + 
        geom_point(aes(x=mean(df$x), y=mean(df$y)), shape=3,colour="red",size=6)
      
      
      if ((pts==TRUE) & (hull==TRUE)) p <- p + geom_point(data=df[df$w==TRUE, ], aes(x, y), shape=21,colour='black') +
        geom_polygon(data = hulls, alpha = 0.1,colour= "gray",fill="gray") 
      
      if ((pts==TRUE) & (hull==FALSE)) p <- p + geom_point(data=df[df$w==TRUE, ], aes(x, y), shape=21,colour='black')
      
      if ((pts==FALSE) & (hull==TRUE)) p <- p + geom_polygon(data = hulls, alpha = 0.1,colour= "gray",fill="gray") 
      
      
      if (plotpdf==TRUE) ggsave(filename = file.path(results_dir,resultdir,paste0("figs"),paste0("pcoa",var,".pdf")),plot=p) else p 
      
    }
    
    a <- pcoa.funk.onevar(data=FR_data,pco=pco_data,var="Din",
                          plotpdf=FALSE,resultdir=taxa,axis.x=2,axis.y=4,
                          jitval=500,pts=TRUE,hull=TRUE,Q1="Q90_D")
    
    b <- pcoa.funk.onevar(data=FR_data,pco=pco_data,var="Rin",
                          plotpdf=FALSE,resultdir=taxa,axis.x=2,axis.y=4,
                          jitval=500,pts=TRUE,hull=TRUE,Q1="Q90_R")
    
    grid.arrange(a,b,ncol=2)
    
  ##Plot pcoa with visualisation of D and hull of DR =  D75R75, D25R25, D25R75, D75R25, D75 or R75
    
    pcoa.funk.dr<-function(data,pco,plotpdf,resultdir,axis.x,axis.y,jitval,var1,var2,Q1,Q2,DR,Funk){
      
      # data=FR_data
      # pco=pco_data
      # resultdir="birds"
      # plotpdf=FALSE
      # axis.x=2
      # axis.y=4
      # jitval=500
      # var1="Din"
      # var2="Rin"
      # Q1="Q75_D"
      # Q2="Q75_R"
      # DR="D75R75"
      # Funk="Din"
      
      df <- data.frame(x = jitter(pco$vectors[,axis.x],jitval),
                       y = jitter(pco$vectors[,axis.y],jitval),
                       z1 = data$FR[rownames(pco$vectors),var1],
                       z2 = data$FR[rownames(pco$vectors),var2])
      
      if (DR=="D75R75") df <- data.frame(df,w=df$z1>as.numeric(data$Q[Q1][[1]]) & df$z2>as.numeric(data$Q[Q2][[1]]))
      if (DR=="D25R25") df <- data.frame(df,w=df$z1<as.numeric(data$Q[Q1][[1]]) & df$z2<as.numeric(data$Q[Q2][[1]]))
      if (DR=="D25R75") df <- data.frame(df,w=df$z1<as.numeric(data$Q[Q1][[1]]) & df$z2>as.numeric(data$Q[Q2][[1]]))
      if (DR=="D75R25") df <- data.frame(df,w=df$z1>as.numeric(data$Q[Q1][[1]]) & df$z2<as.numeric(data$Q[Q2][[1]]))
      
      if (DR=="D75") df <- data.frame(df,w=df$z1>as.numeric(data$Q[Q1][[1]]))
      if (DR=="R75") df <- data.frame(df,w=df$z2>as.numeric(data$Q[Q2][[1]]))
      
      #names(df)[6]<-"w"
      
      df2 <-df[df$w==TRUE,] 
      
      df2 <- df2[complete.cases(df2), ] # needed because there is one NA in the birds dataframe
      
      find_hull <- function(df2) df2[chull(df2$x, df2$y), ]
     
      hulls <- plyr::ddply(df2, "w", find_hull)
      
      p <- ggplot(df, aes(x, y)) +
        geom_point(aes(colour = df$z1))+ scale_colour_gradientn(colours=c("blue","green", "red"),name=Funk) +
        labs(x = paste0("PC",axis.x),y = paste0("PC",axis.y))+ theme_minimal() + ggtitle(DR) +
        geom_point(data=df[df$w==TRUE, ], aes(x, y), shape=21,colour='black') +
        geom_polygon(data = hulls, alpha = 0.1,colour= "black",fill="gray") 
      
      if (plotpdf==TRUE) ggsave(filename = file.path(results_dir,resultdir,paste0("figs"),paste0("pcoa",DR,".pdf")),plot=p) else p 
      
    }
    
####HEAD
    a <- pcoa.funk.dr(data=FR_mammals, pco=pco_mammals, resultdir=taxa,
                      plotpdf=FALSE, axis.x=2, axis.y=3, jitval=500,
                      var1="Din", var2="Rin", Q1="Q75_D", Q2="Q75_R", 
                      DR="D75R75",Funk="Din")
    
    b <- pcoa.funk.dr(data=FR_birds, pco=pco_birds, resultdir=taxa,
                      plotpdf=FALSE, axis.x=1, axis.y=4, jitval=500)
#
    a <- pcoa.funk.dr(data=FR_data, pco=pco_data, resultdir=taxa,
                      plotpdf=FALSE, axis.x=2, axis.y=3, jitval=500,
                      var1="Din", var2="Rin", Q1="Q75_D", Q2="Q75_R", 
                      DR="D75R75",Funk="Din")
    
    b <- pcoa.funk.dr(data=FR_data, pco=pco_data, resultdir=taxa,
                      plotpdf=FALSE, axis.x=2, axis.y=3, jitval=500,
                      var1="Din", var2="Rin", Q1="Q25_D", Q2="Q25_R", 
                      DR="D25R25",Funk="Din")
    grid.arrange(a,b,ncol=2)
#----    
    FR_mammals
    
  # Other version NOT FINISH
    
    load(file=file.path(results_dir,"birds/data_DR_birds.RData"))
      pco<- pco_birds$vectors[rownames(pco_birds$vectors)%in%rownames(data_DR_birds),]
 

 load(file=file.path(results_dir,"mammals/data_DR_mammals.RData"))  
  pco<- pco_mammals$vectors[rownames(pco_mammals$vectors)%in%rownames(data_DR_mammals),]
  
    jitval <- 500
    df <- data.frame(x = jitter(pco[,2],jitval),
                     y = jitter(pco[,3],jitval))
    #df<-merge(df,data_DR_mammals,by="row.names")
     df<-merge(df,data_DR_mammals,by="row.names")
 
    find_hull <- function(df) df[chull(df$x, df$y), ]
    
    df75 <-subset(df,df$DR_class=="D75R75")
    dfAVG<-subset(df,df$DR_class=="AVG")
    df25<-subset(df,df$DR_class=="D25R25")
    
    df75 <- df75[complete.cases(df75), ] # needed because there is one NA in the birds dataframe
    hulls_D75R75 <- find_hull(df75)
    
    dfAVG <- dfAVG[complete.cases(dfAVG), ] # needed because there is one NA in the birds dataframe
    hulls_AVG <- find_hull(dfAVG)
    
    df25 <- df25[complete.cases(df25), ] # needed because there is one NA in the birds dataframe
    hulls_D25R25 <- find_hull(df25)
  
    ggplot(df, aes(x, y)) +
      geom_point(colour = "grey80",shape=21,fill="#4D4D4D1E")+ 
      geom_point(data=df75, colour="orangered",shape=21,fill="#FF45001E")+
      geom_point(data=df25, colour="#E7B800",shape=21,fill="#E7B8001E")+
      geom_point(data=dfAVG, colour="#00AFBB",shape=21,fill="#00AFBB1E")+
      geom_polygon(data = hulls_D75R75, alpha = 0.1,colour= "orangered",fill="white")+
      geom_polygon(data = hulls_AVG, alpha = 0.1,colour= "#00AFBB",fill="white")+
      geom_polygon(data = hulls_D25R25, alpha = 0.1,colour= "#E7B800",fill="white")+
      labs(x = "PC2",y = "PC4")+ theme_minimal()
      
   scale_colour_gradientn(colours=c("blue","green", "red"),name=Funk) +
      labs(x = paste0("PC",axis.x),y = paste0("PC",axis.y))+ theme_minimal() + ggtitle(DR) +
      geom_point(data=df[df$w==TRUE, ], aes(x, y), shape=21,colour='black') +
      geom_polygon(data = hulls, alpha = 0.1,colour= "black",fill="gray")+
      geom_polygon(data = hulls, alpha = 0.1,colour= "black",fill="gray")+
      geom_polygon(data = hulls, alpha = 0.1,colour= "black",fill="gray") 
    

  
  
  data_DR$colsD25R25 <- NA
  data_DR$colsD25R25[data_DR$DR_class=="D25R25"] <- 
    
    data_DR$colsAVG <- NA
    data_DR$colsAVG[data_DR$DR_class=="AVG"] <- 
    
      
      
      
      
      
      load(file=file.path(results_dir,"mammals/data_DR_mammals.RData"))
    
    
    pcoa.funk.dr<-function(data,pco,plotpdf,data_DR,resultdir,axis.x,axis.y,jitval,var1,var2,Q1,Q2,DR,Funk){
      
      # data=FR_mammals
      # pco=pco_mammals
      # resultdir="birds"
      # plotpdf=FALSE
      # axis.x=2
      # axis.y=3
      # jitval=500
      # var1="Din"
      # var2="Rin"
      # Q1="Q75_D"
      # Q2="Q75_R"
      # DR="D75R75"
      # Funk="Din"
      #data_DR =data_DR_mammals
      
      df <- data.frame(x = jitter(pco$vectors[,axis.x],jitval),
                       y = jitter(pco$vectors[,axis.y],jitval),
                       z = data$FR[rownames(pco$vectors),var])
      
      df<-merge(df,data_DR,by="row.names")
      
      find_hull <- function(df) df[chull(df$x, df$y), ]
      
      df75 <-subset(df,df$DR_class=="D75R75")
      dfAVG<-subset(df,df$DR_class=="AVG")
      df25<-subset(df,df$DR_class=="D25R25")
      
      df75 <- df75[complete.cases(df75), ] # needed because there is one NA in the birds dataframe
      hulls_D75R75 <- find_hull(df75)
      
      dfAVG <- dfAVG[complete.cases(dfAVG), ] # needed because there is one NA in the birds dataframe
      hulls_AVG <- find_hull(dfAVG)
      
      df25 <- df25[complete.cases(df25), ] # needed because there is one NA in the birds dataframe
      hulls_D25R25 <- find_hull(df25)
      
      p <- ggplot(df, aes(x, y)) +
        geom_point(aes(colour = df$z))+ scale_colour_gradientn(colours=c("blue","green", "red"),name=Funk) +
        labs(x = paste0("PC",axis.x),y = paste0("PC",axis.y))+ theme_minimal() + ggtitle(taxa) +
        geom_point(data=df[df$w==TRUE, ], aes(x, y), shape=21,colour='black') +
        geom_polygon(data = hulls_D75R75, alpha = 0.1,colour= "orangered", fill="grey",lwd=1.2) +
        geom_polygon(data = hulls_D25R25, alpha = 0.1,colour= "#E7B800", fill="grey", lwd=1.2)+
        geom_polygon(data = hulls_AVG, alpha = 0.1,colour= "#00AFBB", fill="grey", lwd=1.2)
      
      if (plotpdf==TRUE) ggsave(filename = file.path(results_dir,resultdir,paste0("figs"),paste0("pcoa",DR,".pdf")),plot=p) else p 
      
    }
    
    ####HEAD
    a <- pcoa.funk.dr(data=FR_mammals, pco=pco_mammals, resultdir=taxa,data_DR=data_DR_mammals,
                      plotpdf=FALSE, axis.x=2, axis.y=3, jitval=500,
                      var1="Din", var2="Rin", Q1="Q75_D", Q2="Q75_R", 
                      DR="D75R75",Funk="Din")
    b <- pcoa.funk.dr(data=FR_birds, pco=pco_birds, resultdir=taxa,data_DR=data_DR_birds,
                      plotpdf=FALSE, axis.x=2, axis.y=4, jitval=500,
                      var1="Din", var2="Rin", Q1="Q75_D", Q2="Q75_R", 
                      DR="D75R75",Funk="Din")
    