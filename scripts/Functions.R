#FUNCTIONS 

#DUPLICATED2 : modification of duplicated function to extract duplicated lines
duplicated2 <- function(x){ 
  if (sum(dup <- duplicated(x))==0) 
    return(dup) 
  if (class(x) %in% c("data.frame","matrix")) 
    duplicated(rbind(x[dup,],x))[-(1:sum(dup))] 
  else duplicated(c(x[dup],x))[-(1:sum(dup))] 
}


#IDTOTAXO : get the Scientific name from the ID code 

idtotaxo <- function(id,who) #TODO A reprendre !! 
{
  id="sp21260"
  who="mammals"
  
  if (who=="birds") {
    all<- read.csv(file.path(data_dir,"birds/50km/BirdFuncDat.csv"),sep=";")
    all$SpecID <- paste0("ID_",all$SpecID)
    all$Scientific[all$SpecID==id]
  }
  if (who=="mammals") 
    {
    all<- read.csv(file.path(data_dir,"mammals/MamFuncDat.csv"),sep=";") 
    all$MSW3_ID <- paste0("ID_",all$MSW3_ID)
    all$Scientific[all$MSW3_ID==id]
    }
}



#CRAMER : function to compute Cramer correlation between qualitative data 
cramer <- function(x, y) {
  res <- chisq.test(x, y, correct = FALSE)
  chi2 <- as.numeric(res$statistic)
  n <- length(x)
  p <- length(levels(x))
  q <- length(levels(y))
  m <- min(p - 1, q - 1)
  V <- sqrt(chi2/(n * m))
  return(V)
}


#DISTINCTIVNESS_GLOB : compute the distinctivness of the entire species pool 

  distinctiveness_glob <- function(com_dist,abund){
  
  com_df <- as.data.frame(colnames(com_dist))
  if (!is.null(dim(com_dist))) {
    if (is.null(abund)) {
      num = colSums(com_dist)
      denom = nrow(com_df) - 1
    }
    else {
      num = apply(com_dist, 2, function(x) sum(x * com_df[, 
                                                          abund]))
      denom = sum(com_df[[abund]]) - com_df[[abund]]
    }
  } else {
    denom = 0
  }
  if (length(denom) > 1) {
    com_df[, "Di"] = as.numeric(num/denom)
  } else if (length(denom) == 1 & denom != 0) {
    com_df[, "Di"] = as.numeric(num/denom)
  } else {
    com_df[, "Di"] = NaN
  }
  return(com_df)
}


#ELLIP corelation with ellipse visualisation 

  ellip <- function(data,varsub,plotpdf,ordvar,taxa){
  #varsub <- c("TD_sp","FEve","FDis","FOri","FDiv","FRic","FSpe","UiMean","UiSkw","DiMean","DiSkw","RiMean","RiSkw","Ui90","RFR_Ui90","Ui90_perc")
  #final_results <- final_results_mamals
  #plotpdf=TRUE
  
  require(RColorBrewer)
  require(ellipse)
  
  data_cor=cor(data[,varsub],use="pairwise.complete.obs")
  my_colors <- brewer.pal(5, "Spectral")
  my_colors=colorRampPalette(my_colors)(100)
  ord <- order(data_cor[ordvar, ])
  data_ord <-  data_cor[ord, ord]
  data_ord <- data_cor
  
  if (plotpdf==TRUE) {
    pdf(file.path(results_dir,"mamals/funk_mamals.cor.pdf"))
    plotcorr(data_ord , col=my_colors[data_ord*50+50] , mar=c(1,1,1,1),main=taxa)
    dev.off()
  } else {plotcorr(data_ord , col=my_colors[data_ord*50+50] , mar=c(1,1,1,1),main=taxa)}
  }
  

#G.PLOT : compute ggplot between x and y 
  g.plot <- function(dataplot,x,y,taxa) {
  ggplot(dataplot, aes(get(x), get(y))) + geom_point() + geom_smooth(method = "lm")+ labs(x = x,y=y)  + ggtitle(taxa)
}

#MAP.FUNK : Draw map for each indices for each taxa 
map.Funk <- function(data,map,var,nlevels,plotpdf,resultdir,dalto){
    
    require(viridis)
    require(RColorBrewer)
    
    # data<- funk_mammals
    # var <- "D75R75"
    #nlevels <- 6 #Choix du nombre de classe
    #map<- map_mamals
    #resultdir="mamals"
    # plotpdf=TRUE
    # dalto: true of false
    
    data<- data[,var]
    data[data==0]<-NA
    data<-as.numeric(data)
    if(dalto==TRUE) colour <- viridis(nlevels) else{ colour <- rev(brewer.pal(nlevels, "Spectral"))}

    bInf <- range(na.omit(data))[1]
    bSup <- range(na.omit(data))[2]
    vect = seq(bInf,bSup,length.out=(nlevels+1))
    RSLAB = cut(data, breaks=vect,include.lowest=TRUE,right=FALSE)
    cols <- colour[RSLAB]
    cols[is.na(cols)] <- "white" #all NAs will be plot in black
    colsborder<-cols
    colsborder[colsborder=="white"]<-"gray82"
    
    if (plotpdf==TRUE){
      pdf(file.path(results_dir,resultdir,paste0("figs"),paste0("map",var,".pdf")))
      nf=layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
      plot(map,col=cols,main=var,lwd=0.005,border=colsborder)
      hist(data,main=paste0("Histogram of ",var),xlab="", las=2)
      barplot(summary(RSLAB)[-(nlevels+1)],col=c(colour),ylab="Number of cells", las=2, cex.names = 0.6)
      dev.off()
    } else {  
      nf=layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
      plot(map,col=cols,main=var,lwd=0.005,border=colsborder)
      hist(data,main=paste0("Histogram of ",var),xlab="", las=2)
      barplot(summary(RSLAB)[-(nlevels+1)],col=c(colour),ylab="Number of cells", las=2, cex.names = 0.6)}
      }
  
#PANEL COR : Function for correlation graphics 
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
  {
    usr <- par("usr"); on.exit(par(usr)) 
    par(usr = c(0, 1, 0, 1)) 
    r <- abs(cor(x, y)) 
    txt <- format(c(r, 0.123456789), digits=digits)[1] 
    txt <- paste(prefix, txt, sep="") 
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
    
    test <- cor.test(x,y) 
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " ")) 
    
    text(0.5, 0.5, txt, cex = cex * r) 
    text(.8, .8, Signif, cex=cex, col=2) 
  }  
  
#MAP.FUNK2 nice fig 
  
  map.Funk <- function(data,map,var,nlevels,plotpdf,resultdir,dalto){
    
    require(viridis)
    require(RColorBrewer)
    
    # data<- funk_mammals
    # var <- "D75R75"
    #nlevels <- 6 #Choix du nombre de classe
    #map<- map_mamals
    #resultdir="mamals"
    # plotpdf=TRUE
    # dalto: true of false
    
    data<- data[,var]
    data[data==0]<-NA
    data<-as.numeric(data)
    if(dalto==TRUE) colour <- viridis(nlevels) else{ colour <- rev(brewer.pal(nlevels, "Spectral"))}
    
    bInf <- range(na.omit(data))[1]
    bSup <- range(na.omit(data))[2]
    vect = seq(bInf,bSup,length.out=(nlevels+1))
    RSLAB = cut(data, breaks=vect,include.lowest=TRUE,right=FALSE)
    cols <- colour[RSLAB]
    cols[is.na(cols)] <- "white" #all NAs will be plot in black
    colsborder<-cols
    colsborder[colsborder=="white"]<-"gray82"
    
    
    if (plotpdf==TRUE){
      pdf(file.path(results_dir,resultdir,paste0("figs"),paste0("map",var,".pdf")))
      nf=layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
      plot(map,col=cols,main=var,lwd=0.005,border=colsborder)
      
      
      
  library(viridis)
  p <- ggplot() +
    geom_polygon(data = spdf_fortified, aes(fill = nb_equip, x = long, y = lat, group = group) , size=0, alpha=0.9) +
    theme_void() +
    scale_fill_viridis(trans = "log", breaks=c(1,5,10,20,50,100), name="Number of restaurant", guide = guide_legend( keyheight = unit(3, units = "mm"), keywidth=unit(12, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) ) +
    labs(
      title = "South of France Restaurant concentration",
      subtitle = "Number of restaurant per city district", 
      caption = "Data: INSEE | Creation: Yan Holtz | r-graph-gallery.com"
    ) +
    theme(
      text = element_text(color = "#22211d"), 
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      panel.background = element_rect(fill = "#f5f5f2", color = NA), 
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      
      plot.title = element_text(size= 22, hjust=0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
      plot.subtitle = element_text(size= 17, hjust=0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
      plot.caption = element_text( size=12, color = "#4e4d47", margin = margin(b = 0.3, r=-99, unit = "cm") ),
      
      legend.position = c(0.7, 0.09)
    ) +
    coord_map()
  p 
    }} # TODO j'ai ajouté ca mais est ce bon ? 
    
    
#CLEAN NAME  
  stringCleaning <- function(x) {
    #   x <- stringr::str_trim(x)
    #   x <- tolower(x)
    #   x <- gsub("\\s+", " ", x)
    #   x <- gsub("[^[:space:]A-z0-9]", "", x)
    stringr::str_trim(tolower(gsub("\\s+", " ", gsub("[^[:space:]A-z0-9]", "", x))))
}
  
  
#TEST.COREL : Function to test correlation between traits
  # cramer between qualitative vars
  # Pearson between quantitative vars
  # lm between quantitative & qualitative vars
  test.corel <- function(traits){
    res.test.corel <- list() 
    
    quanqual<- sapply(traits,is.numeric)
    quanti<-names(quanqual)[quanqual%in%TRUE]
    quali<-names(quanqual)[quanqual%in%FALSE]
    
    #Entre traits quanti
    quanti_traits<-data.frame(traits[,quanti])
    res.test.corel[[1]]<-cor(quanti_traits)
    
    #Entre traits quali
    combinaison<- combn(quali,m=2)
    res.test.corel[[2]] <- do.call(rbind,lapply(1:dim(combinaison)[2],function(i){
      r=cramer(traits[,combinaison[1,i]],traits[,combinaison[2,i]])
      cbind.data.frame(t1=combinaison[1,i],t2=combinaison[2,i],r=r)
    }))
    
    #Entre traits quali et quanti
    corqalquan <- do.call(rbind,lapply(quali,function(IDqual){
      #IDqual <- quali[1]
      cor1 <- do.call(cbind,lapply(quanti,function(IDquan){
        summary(lm(as.formula(paste(IDquan,"~",IDqual)),data=traits[,c(IDqual,IDquan)]))[[8]]
      }))
    })) 
    colnames(corqalquan) <- quanti
    rownames(corqalquan) <- quali
    res.test.corel[[3]] <- corqalquan
    
    names(res.test.corel) <- c("quanquan","qualqual","quanqual")
    return(res.test.corel)
  }  
  
  
#TODO : commentaires 
  who.remote <- function(who,remote){
  if (remote %in% TRUE) {
    data_dir <<- file.path("/home","nmouquet","Dropbox","SCIENCE","RALL","RALLL_R","data")
    results_dir <<- file.path("/home","nmouquet","Dropbox","SCIENCE","RALL","RALLL_R","results")
    script_dir <<- file.path("/home","nmouquet","Dropbox","SCIENCE","RALL","RALLL_R","scripts")
  } else {
    if (who %in% "NM") {
      data_dir <<- "data"
      results_dir <<- "results"
    } else {
      data_dir <<- file.path("~/Dropbox","RALL","RALLL_R","data")
      results_dir <<- file.path("~/Dropbox","RALL","RALLL_R","results")
    }
    script_dir <<-"scripts"
  }
}

