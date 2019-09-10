library(raster)
root      <- "~/Dropbox/RALL/RALLL_R/FUNCRARITY/"
src_mammals <- subset(src, src$taxa.x=="mammals")
src_mammals_D75R75 <- subset(src_mammals,src_mammals$DR_class=="D75R75")

species = "sp41026"

spname = subset(taxaInfo_mammals,taxaInfo_mammals$ID == species)


par(mfrow=c(1,2))
cur = get(load(paste0("/Volumes/ROBERTO/ahasverus/outputs/projs/proj_", species, "_1979-2013_bins")))
fut = get(load(paste0("/Volumes/ROBERTO/ahasverus/outputs/projs/proj_", species, "_2061-2080_bins")))
plot(cur)
title(spname$Name)
plot(fut)

fls = list.files("/Volumes/ROBERTO/", recursive = TRUE, pattern = "_1979-2013_bins$", full.names = TRUE)
fls = list.files("/Volumes/ROBERTO/", recursive = TRUE, pattern = "^proj.+_2061-2080_bins$", full.names = TRUE)

for (i in 1:length(fls)) {
  cat(i, "\r")
  tmp = get(load(fls[i]))
  if (i == 1) {
    rich = tmp
  }else{
    rich[] = rich[] + tmp[]
  }
}

plot(rich)
