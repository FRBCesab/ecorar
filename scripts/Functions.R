#FUNCTIONS 
#----
# COLOR TERMINAL BRANCHES: function to color only the last branch from Nick Crouch's GITHUB
# phy: object of class phylo
# data: named numeric vector
# breaks: number of color breaks
# cols: character vector of length 2, start col, end col
# non.terminal.col: vector character length 1, color non terminal branches

color.terminal.branches <- function(phy, data, breaks=6, cols=c("black", "red"), non.terminal.col= "black", edge.width=1, col.bias=1, legend.title = "", show.tip.label=FALSE, alt.col.data=NULL,...){
  
  sapply(phy$tip.label, check.spp.search.duplicates, phy$tip.label)
  
  # Number of tips in the phylogeny
  n.tips <- length(phy$tip.label)
  
  # This will be the main vector which colors will be assigned to
  col.vector <- vector(mode="character",length=nrow(phy$edge))
  
  # Non-terminal branches will be black
  col.vector[phy$edge[,2]>n.tips] <- non.terminal.col
  
  
  # Create a vector of colors 
  colors <-  colorRampPalette(cols, bias=col.bias)(breaks)
  
  # range of trait values present
  if(is.null(alt.col.data)==TRUE){
    range.vals <- range(data)
  } else {
    range.vals <- range(c(alt.col.data, data))
  }
  diff <- range.vals[2] - range.vals[1]
  
  bin.size <- diff / breaks
  
  bin.ref <- as.data.frame(matrix(NA, ncol=3, nrow=breaks))
  colnames(bin.ref) <- c("Bin","Start","End")
  bin.ref$Bin <- seq(1, breaks, 1)
  bin.ref[1,2] <- range.vals[1]
  bin.ref[1,3] <- range.vals[1] + bin.size
  
  for(i in 2:nrow(bin.ref)){
    
    bin.ref[i,2] <- bin.ref[i-1,3]+0.00001
    bin.ref[i,3] <- bin.ref[i,2] + bin.size
    
  }
  
  
  ## for legend
  legend.text <- vector()
  
  for(i in 1:nrow(bin.ref)){
    
    t <- paste(as.character(bin.ref[i,2]), " - ",as.character(bin.ref[i,3]))
    
    legend.text[i] <- t
    
  }
  
  
  # create a duplicate of phy$edge which can be manipulated as required
  edge.data <- as.data.frame(phy$edge)
  
  
  # Yes, its a loop, but its easier
  for(i in 1:length(phy$tip.label)){
    
    spp <- phy$tip.label[i]
    
    spp.trait.val <- as.numeric(data[grep(spp,names(data))])
    
    bin.num <- bin.ref[spp.trait.val >= bin.ref$Start & spp.trait.val <= bin.ref$End,1]
    
    spp.col <- colors[bin.num]
    
    edge.row <- as.numeric(rownames(edge.data[edge.data$V2==i,]))
    
    col.vector[edge.row] <- spp.col
    
  }
  
  max.val <- range.vals[2]
  
  #layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
  #l <-  matrix(c(1,1,2),ncol=3,nrow=1)
  #layout(l)
  #par(mfrow=c(1,2))
  #par(mai=rep(0.1,4))
  plot(phy, show.tip.label=show.tip.label, edge.col=col.vector, edge.width=edge.width,type="fan",tip.color="white")
  
  #image <- matrix(colors, ncol=1)
  #image <- as.matrix(rev(legend_image))
  #legend_image <- as.raster(image)
  # legend_image <- as.raster(matrix(rev(colors), ncol=1))
  
  # plot(c(0,.9),c(0,max.val),type = 'n', axes = F,xlab = '', ylab = '', main = legend.title, xlim=c(0,.2))
  
  # rasterImage(legend_image, 0, 0, .08,max.val)
  # text(x=0.11, y = seq(0,max.val,l=5), labels = seq(0,max.val,l=5))
}

###

check.spp.search.duplicates <- function(search.spp, all.spp){
  
  returns <- grep(search.spp, all.spp)
  
  if(length(returns)>1){
    
    stop("The species  ",search.spp,"  will cause problems - searching for this species in all the species present will return multiple species, consider renaming!")
    
  }
  
}

#----

# COLOR TERMINAL BRANCHES: function to compute only the last branch from Nick Crouch's GITHUB
# phy:        Phylogeny of class 'phylo'
# data:       'data.frame' containing data, must have a 'Species' column
# plot.data:  column in data to be plotted. Must be numeric, starting at 1
# pch:        plotting symbol
# cols:       vector of cols?

add.tip.icons <- function(phy, data, plot.data="Location", grouping="Species", pch=21, cols=c("red","blue","green"), cex=2){
  
  require(ape)
  
  if( (plot.data %in% colnames(data)) == FALSE){
    stop("Specified data to plot does not appear in the column names of the data")
  }
  
  v <- data[,grep(grouping, colnames(data))] %in% phy$tip.label
  
  if((FALSE %in% v) == TRUE){
    stop("Tip labels of phylogeny do not match grouping vector specified")
  }
  
  
  d <- data[,grep(plot.data, colnames(data))]
  
  # does number of unique entries in plot.data equal number of colors provided?
  num.unique.vals <- length(unique(d))
  
  if(num.unique.vals != length(cols)){
    
    stop("Number of states does not equal number of colors provided")
    
  }
  
  names(d) <- data[,grep(grouping, colnames(data))]
  
  names <- phy$tip.label
  
  m <- match(names, names(d))
  
  r <- d[c(m)]
  
  plot.col <- vector(length=length(phy$tip.label))
  
  for(i in 1:length(cols)){
    
    plot.col[as.numeric(r)==i] <- cols[i]
    
  }
  
  
  # Plot
  
  plot(phy, show.tip.label=FALSE)
  
  if(pch < 21){
    
    tiplabels(col=plot.col, pch=pch, cex=cex)
    
  } else if(pch >= 21 & pch <= 25){
    
    tiplabels(bg=plot.col, pch=pch, cex=cex)
    
  } else {
    
    stop("Invalid plotting symbol value provided")
    
  }
  
  
}
#----

#ARC.CLADELABEL : add arc around phylogenie tree
#Modifer version of ARC.CLADELABEL from phytools packages that do not allow change in color of arc 
arc.cladelabels<-function(tree=NULL,text,node,ln.offset=1.02,
                          lab.offset=1.06,cex=1,orientation="curved",colarc="black",coltext="black",...){
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(obj$type!="fan") stop("method works only for type=\"fan\"")
  h<-max(sqrt(obj$xx^2+obj$yy^2))
  if(hasArg(mark.node)) mark.node<-list(...)$mark.node
  else mark.node<-TRUE
  if(mark.node) points(obj$xx[node],obj$yy[node],pch=21,
                       bg="red")
  if(is.null(tree)){
    tree<-list(edge=obj$edge,tip.label=1:obj$Ntip,
               Nnode=obj$Nnode)
    class(tree)<-"phylo"
  }
  d<-phytools::getDescendants(tree,node)
  d<-sort(d[d<=Ntip(tree)])
  deg<-atan(obj$yy[d]/obj$xx[d])*180/pi
  ii<-intersect(which(obj$yy[d]>=0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]>=0))
  deg[ii]<-360+deg[ii]
  plotrix::draw.arc(x=0,y=0,radius=ln.offset*h,deg1=min(deg),
           deg2=max(deg),col=colarc, ...)
  if(orientation=="curved")
    plotrix::arctext(text,radius=lab.offset*h,
            middle=mean(range(deg*pi/180)),cex=cex,col=coltext,...)
  else if(orientation=="horizontal"){
    x0<-lab.offset*cos(median(deg)*pi/180)*h
    y0<-lab.offset*sin(median(deg)*pi/180)*h
    text(x=x0,y=y0,label=text,
         adj=c(if(x0>=0) 0 else 1,if(y0>=0) 0 else 1),
         offset=0)
  }
}
  
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
    #resultdir="mammals"
    # plotpdf=TRUE
    # dalto: true of false
  #  log10
  
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
      barplot(summary(RSLAB)[-(nlevels+1)],col=c(colour),ylab="Number of cells", las=2, cex.names = 0.6)
      }
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

  
  
  
#TODO : TWO FUNCTION TO CHANGE COLOR OF FIGURE IN PHYLOGENY PLOT   
  # function to calculate brightness values
  brightness <- function(hex) {
    v <- col2rgb(hex)
    sqrt(0.299 * v[1]^2 + 0.587 * v[2]^2 + 0.114 * v[3]^2) /255
  }
  # given a color ramp, map brightness to ramp also taking into account 
  # the alpha level. The defaul color ramp is grey
  #
  img_to_colorramp <- function(img, ramp=grey) {
    cv <- as.vector(img)
    b <- sapply(cv, brightness)
    g <- ramp(b)
    a <- substr(cv, 8,9)     # get alpha values
    ga <- paste0(g, a)       # add alpha values to new colors
    img.grey <- matrix(ga, nrow(img), ncol(img), byrow=TRUE)  
  }
