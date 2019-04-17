load(file=file.path(results_dir,"birds","50km","Di_locall_birds.RData"))
mean_loc_Di <- apply(Di_locall,1,function(x) mean(x, na.rm = T)) 
mean_loc_Di <- data.frame(mean_loc_Di)

mean_loc_Di <- apply(Di_locall,2,function(x) mean(x, na.rm = T)) 
SP<-funk_birds[funk_birds$cell %in% names(mean_loc_Di),]$TD_sp
plot(mean_loc_Di~SP)


Di_loc_glob <- merge(Di,mean_loc_Di,by="row.names",all.x=FALSE)
colnames(Di_loc_glob) <- c("ID","Di_glob","Di_loc")

ggplot(Di_loc_glob, aes(x=Di_glob, y=Di_loc)) + 
  geom_point(size=1) + 
  stat_smooth(method = "lm", formula = y ~ x, size = 1,se=TRUE) + 
  ggtitle("Birds") + 
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1)

#Extract outliers species 
Di_loc_glob <- Di_loc_glob[Di_loc_glob$ID %in% birdsID$ID,]


up <- Di_loc_glob[Di_loc_glob[,"Di_glob"]<0.55,]
up <- up[up[,"Di_loc"]>0.7,]
upspecies <- birdsID[birdsID$ID %in% up$ID,]

down <- Di_loc_glob[Di_loc_glob[,"Di_glob"]>0.75,]
down <- down[down[,"Di_loc"]<0.5,]
downspecies <- birdsID[birdsID$ID %in% down$ID,]


ids <- upspecies$ID
subUP <- mclapply(ids,function(id) {  
  spe_up <-  names(unlist(lapply(occ_birds_list, function(x) x[x==id])))
},mc.cores = 4)
names(subUP) <- ids
subUP <- data.frame(ID=unique(unlist(subUP)),UP=rep(1, length(unique(unlist(subUP)))))


ids <- downspecies$ID
subdown <- mclapply(ids,function(id) {  
  spe_down <-  names(unlist(lapply(occ_birds_list, function(x) x[x==id])))
},mc.cores = 4)
names(subdown) <- ids
subdown <- data.frame(ID=unique(unlist(subdown)),DOWN=rep(1, length(unique(unlist(subdown)))))


RSsubdown <- funk_birds[funk_birds$cell %in% subdown$ID,]
RSsubUP <- funk_birds[funk_birds$cell %in% subUP$ID,]

boxplot(RSsubdown$TD_sp,RSsubUP$TD_sp,funk_birds$TD_sp,col=c("green","red","grey"))




map<-readOGR(file.path(data_dir,"ReferenceGrid50Km","gridLand50km.shp"))
#names of each cell
ID_cell<-map@data[,1]


map <- merge(map, subUP, by = 'ID')
map <- merge(map, subdown, by = 'ID')
map@data[is.na(map@data )]<-0

pdf(file.path(results_dir,"birds","50km",paste0("figs"),paste0("map","testDILocUP_Map",".pdf")))
spplot(map["UP"],col.regions = c("gray88","red"),
              ## set the border color and width
              col="transparent",
              par.settings = list(axis.line=list(col="black")),
              contour = T) 
dev.off()



pdf(file.path(results_dir,"birds","50km",paste0("figs"),paste0("map","testDILocDOWN_Map",".pdf")))
spplot(map["DOWN"],col.regions = c("gray88","green"),
       ## set the border color and width
       col="transparent",
       par.settings = list(axis.line=list(col="black")),
       contour = T) 

dev.off()

