

load(file=file.path(results_dir,"mammals/mammalsID.RData"))
load(file=file.path(results_dir,"mammals/mammalstrait.RData"))

mammalsID$clean <- paste(word(mammalsID$checkname,1),word(mammalsID$checkname,2),sep=" ")
mammalstrait$clean <- paste(word(mammalstrait$checkname,1),word(mammalstrait$checkname,2),sep=" ")

mammalsID$clean <- stringCleaning(mammalsID$clean)
mammalstrait$clean <- stringCleaning(mammalstrait$clean)


dim(mammalsID[mammalsID$clean %in% mammalstrait$clean,])
dim(mammalstrait[mammalstrait$clean %in% mammalsID$clean,])


df1 <- tbl_df(mammalsID)
df2 <- tbl_df(mammalstrait)


df3<- df1 %>%
  stringdist_inner_join(df2, by =  "clean",max_dist=2)

as.data.frame(df3)



#---
##############SUR LES MAT DE DEPARt
duplicated2 <- function(x){ 
  if (sum(dup <- duplicated(x))==0) 
    return(dup) 
  if (class(x) %in% c("data.frame","matrix")) 
    duplicated(rbind(x[dup,],x))[-(1:sum(dup))] 
  else duplicated(c(x[dup],x))[-(1:sum(dup))] 
}
############



load(file=file.path(data_dir,"mammals/MammalsID.RData"))
mammalsID<-data.frame(Mammals)

mammalstrait_all<- read.csv(file.path(data_dir,"mammals","MamFuncDat.csv"),sep=";") 

mammalstrait_all$Diet.Inv <- as.numeric(mammalstrait_all$Diet.Inv)
mammalstrait_all$Diet.Vend <- as.numeric(mammalstrait_all$Diet.Vend)
mammalstrait_all$Diet.Vect <- as.numeric(mammalstrait_all$Diet.Vect)
mammalstrait_all$Diet.Vfish <- as.numeric(mammalstrait_all$Diet.Vfish)
mammalstrait_all$Diet.Vunk <- as.numeric(mammalstrait_all$Diet.Vunk)
mammalstrait_all$Diet.Scav <- as.numeric(mammalstrait_all$Diet.Scav)
mammalstrait_all$Diet.Fruit <- as.numeric(mammalstrait_all$Diet.Fruit)
mammalstrait_all$Diet.Nect <- as.numeric(mammalstrait_all$Diet.Nect)
mammalstrait_all$Diet.Seed <- as.numeric(mammalstrait_all$Diet.Seed)
mammalstrait_all$Diet.PlantO <- as.numeric(mammalstrait_all$Diet.PlantO)
mammalstrait_all$ForStrat.Value<- as.factor(mammalstrait_all$ForStrat.Value)
mammalstrait_all$Activity.Nocturnal<- as.factor(mammalstrait_all$Activity.Nocturnal)
mammalstrait_all$Activity.Crepuscular<- as.factor(mammalstrait_all$Activity.Crepuscular)
mammalstrait_all$Activity.Diurnal<- as.factor(mammalstrait_all$Activity.Diurnal)
mammalstrait_all$BodyMass.Value<- as.numeric(mammalstrait_all$BodyMass.Value)

interest.var.mammals <- c("MSW3_ID","Scientific","Diet.Inv","Diet.Vend","Diet.Vect","Diet.Vfish","Diet.Vunk","Diet.Scav","Diet.Fruit","Diet.Nect","Diet.Seed",
                          "Diet.PlantO","ForStrat.Value","Activity.Nocturnal","Activity.Crepuscular","Activity.Diurnal","BodyMass.Value")

mammalstrait <- mammalstrait_all[,interest.var.mammals]


mammalsID$clean <- paste(word(mammalsID$Name,1),word(mammalsID$Name,2),sep=" ")
mammalstrait$clean <- paste(word(mammalstrait$Scientific,1),word(mammalstrait$Scientific,2),sep=" ")

#clean name
mammalsID$clean <- stringCleaning(mammalsID$clean)
mammalstrait$clean <- stringCleaning(mammalstrait$clean)

#8 species are dupblicated due to subspecies (no functional trait for them, delete)

#Find Species
sp_dup <- mammalsID[duplicated2(mammalsID[, "clean"]), ]
#Choose to delete spp
sp_dup<-sp_dup[grep("ssp.",sp_dup$Name, ignore.case = TRUE),]
mammalsID <- mammalsID[!rownames(mammalsID)%in%rownames(sp_dup),]

dim(mammalsID[mammalsID$clean %in% mammalstrait$clean,])
dim(mammalstrait[mammalstrait$clean %in% mammalsID$clean,])


mammalstrait$checkname<-NA
for (i in 1:length(mammalstrait$clean)){
  mammalstrait$checkname[i] <- as.matrix(gnr_resolve(names = as.character(mammalstrait$clean[i])))[1,3]
  print(paste0("i",i))
}


mammalsID$checkname<-NA
for (i in 1:length(mammalsID$Name)){
mammalsID$checkname[i] <- as.matrix(gnr_resolve(names = as.character(mammalsID$Name[i])))[1,3]
print(paste0("i",i))
}
#mammalsID$checkname[5128] <- as.character("Pseudoberylmys muongbangensis")

#Converte occ_mammals in list
occ_mammals_list <- lapply(occ_mammals, function(.col){
  sp <- row.names(occ_mammals)[.col == 1] })





df1 <- tbl_df(mammalsID)
df2 <- tbl_df(mammalstrait)

df3<- df1 %>%
  stringdist_inner_join(df2, by =  "clean",max_dist=2)

dim(as.data.frame(df3))

