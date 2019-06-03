library(ggplot2)
library(grid)
library(gridExtra)

ras
library(raster)
library(tiff)

footprint=raster(file.path(data_dir,"GlobalHumanFootprintIndex.tiff"))
ras

footprint <- raster(footprint)
new.rasterA = projectRaster(footprint, ras) #define the projection and extent

r.stack = stack(new.rasterA, rasterB)

new_2013 <- crop(extend(footprint, ras), ras)
all.equal(extent(raster_2015), extent(new_2013))
pco2<- vegan::wcmdscale(disTraits_birds) 



load(file=file.path(results_dir,"birds/pco_birds.RData"))
load(file=file.path(results_dir,"mammals/pco_mammals.RData"))

#Function to add levels
addLevel <- function(x, newlevel=NULL) {
  if(is.factor(x)) {
    if (is.na(match(newlevel, levels(x))))
      return(factor(x, levels=c(levels(x), newlevel)))
  }
  return(x)
}



pcoa.funk.dr<-function(data,pco,plotpdf,resultdir,axis.x,axis.y,jitval,var1,var2,Q1,Q2,DR,Funk){
  
  pco <- pco_mammals
  traits <- mammalstrait
  data <- FR_mammals
  taxa <- "mammals"
  DR <-data_DR_mammals
  
  pco <- pco_birds
  pco <- pco2
  traits <- birdstrait
  data <- FR_birds
  taxa <- "birds"
  DR <-data_DR_birds
  
  var1="Din"
  var2="Rin"
  axis.x=1
  axis.y=2
  jitval=500
  df <- data.frame(x = jitter(pco$vectors[,axis.x],jitval),
                   y = jitter(pco$vectors[,axis.y],jitval),
                   z1 = data$FR[rownames(pco$vectors),var1],
                   z2 = data$FR[rownames(pco$vectors),var2])
  df <- merge(df,DR,by ="row.names")
  
  df <- data.frame(x = jitter(pco2[,axis.x],jitval),
                   y = jitter(pco2[,axis.y],jitval),
                   z1 = data$FR[rownames(pco2),var1],
                   z2 = data$FR[rownames(pco2),var2])
  df <- merge(df,DR,by ="row.names")
  
  
  dfD75R75 <- subset(df,df$DR_class=="D75R75")
  dfD25R25 <- subset(df,df$DR_class=="D25R25")
  dfAVG <- subset(df,df$DR_class=="AVG")

  #names(df)[6]<-"w"
  dfD75R75_2 <- dfD75R75[complete.cases(dfD75R75), ] # needed because there is one NA in the birds dataframe
  find_hull_D75R75 <- function(dfD75R75_2) dfD75R75_2[chull(dfD75R75_2$x, dfD75R75_2$y), ]
  hulls_D75R75  <- plyr::ddply(dfD75R75_2, "DR_class", find_hull_D75R75)
  
  dfD25R25_2 <- dfD25R25[complete.cases(dfD25R25), ] # needed because there is one NA in the birds dataframe
  find_hull_D25R25 <- function(dfD25R25_2) dfD25R25_2[chull(dfD25R25_2$x, dfD25R25_2$y), ]
  hulls_D25R25  <- plyr::ddply(dfD25R25_2, "DR_class", find_hull_D25R25)
  
  dfAVG_2 <- dfAVG[complete.cases(dfAVG), ] # needed because there is one NA in the birds dataframe
  find_hull_AVG <- function(dfAVG_2) dfAVG_2[chull(dfAVG_2$x, dfAVG_2$y), ]
  hulls_AVG  <- plyr::ddply(dfAVG_2, "DR_class", find_hull_AVG)
  

  df$DR_class <- addLevel(df$DR_class, "Others")
  df$DR_class[is.na(df$DR_class)] <-"Others"
  
  cols <- c("AVG" = "#00AFBB", "D25R25" = "#E7B800", "D25R75" = "#CCCCCC66", "D75R25" = "#CCCCCC66", "D75R75"= "orangered","Others" = "#CCCCCC66")

 a <- ggplot(df, aes(x, y)) + ylim(-0.7,0.4) + xlim(-0.5,0.5) +
    geom_point(fill="#E0E0E033", shape=21, size = 3, stroke = 1, aes(colour = df$DR_class)) + 
    labs(x = paste0("PC",axis.x),y = paste0("PC",axis.y))+ theme_minimal() + ggtitle(taxa) + 
    scale_colour_manual(values=cols) + theme_bw()+scale_fill_manual(values=cols) + 
    geom_polygon(data = hulls_D75R75, alpha = 0.1,colour= "orangered",fill="white",lwd=1) + 
    geom_polygon(data = hulls_D25R25, alpha = 0.1,colour= "#E7B800",fill="white",lty=2,lwd=1) +
    geom_polygon(data = hulls_AVG, alpha = 0.1,colour= "#00AFBB",fill="white",lty=3,lwd=1) +
    theme(legend.position="none")
    
 b <- ggplot(df, aes(x, y)) + ylim(-0.7,0.4) + xlim(-0.5,0.5) +
    geom_point( shape=21, size = 3, stroke = 1, aes(colour = df$DR_class,fill=df$DR_class,alpha=0.4)) + 
    labs(x = paste0("PC",axis.x), y = paste0("PC",axis.y))+ theme_minimal() + ggtitle(taxa) + 
    scale_colour_manual(values=cols) + theme_bw()+scale_fill_manual(values=cols) + 
    geom_polygon(data = hulls_D75R75, alpha = 0.1,colour= "orangered",fill="white",lwd=1) + 
    geom_polygon(data = hulls_D25R25, alpha = 0.1,colour= "#E7B800",fill="white",lty=2,lwd=1) +
    geom_polygon(data = hulls_AVG, alpha = 0.1,colour= "#00AFBB",fill="white",lty=3,lwd=1) +
    theme(legend.position="none")
    
 grid.arrange(a,b,ncol=2)   
}


makeTransparent("grey80",alpha=0.4)

data_DR$colsD25R25 <- NA
data_DR$colsD25R25[data_DR$DR_class=="D25R25"] <- "#"

data_DR$colsAVG <- NA
data_DR$colsAVG[data_DR$DR_class=="AVG"] <- 
