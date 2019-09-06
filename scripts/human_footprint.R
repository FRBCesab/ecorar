library(sp)
library(rgdal)
library(rgeos)
library(raster)

path_data <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/data/"           ### !!!

study <- raster(paste0(path_data, "reference_grid_50km.tif"))

foot  <- raster(paste0(path_data, "HFP2009.tif"))
foot  <- projectRaster(from = foot, to = study)

coords <- xyFromCell(foot, 1:length(foot))
values <- foot[]

dat <- data.frame(coords, values)
dat <- dat[!is.na(dat$value), ]

dat <- SpatialPointsDataFrame(coords = dat[ , 1:2], data = data.frame(value = dat[ , 3]))
proj4string(dat) <- proj4string(study)

hfp <- rasterize(x = dat, y = study, field = "value", fun = mean)

writeRaster(x = hfp, filename = "HFP_2009.tif", format = "GTiff", overwrite = TRUE)
