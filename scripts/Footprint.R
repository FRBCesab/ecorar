rm(list=ls(all=TRUE)) 
source("./scripts/Functions.R")
who.remote(remote=FALSE,who="NL")
library(rgeos)
library(ggplot2)
library(grid)
library(gridExtra)
library(raster)
library(tiff)

<<<<<<< HEAD
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
=======
footprint=raster(file.path(data_dir,"GlobalHumanFootprintIndex.tiff"))

>>>>>>> 1f7759196551b79bfbf109d5fff121ac449b6a78

footprint <- raster(footprint)
crs(footprint)="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs "

map_shp<-readOGR(file.path(data_dir,"ReferenceGrid50Km","gridLand50km.shp")) 
map_extent      <- map_shp@bbox
poly_coords     <- lapply(map_shp@polygons, function(x) x@Polygons[[1]]@coords)
grid_res_x      <- unique(unlist(lapply(poly_coords, function(x) max(x[ , 1]) - min(x[ , 1]))))
grid_res_y      <- unique(unlist(lapply(poly_coords, function(x) max(x[ , 2]) - min(x[ , 2]))))

cells_infos <- data.frame(
  cell_id = map_shp@data[ , 1],
  gCentroid(map_shp, byid = TRUE)
)

ras <- raster(
  xmn         = map_extent["x", "min"],
  xmx         = map_extent["x", "max"],
  ymn         = map_extent["y", "min"],
  ymx         = map_extent["y", "max"],
  crs         = proj4string(map_shp),
  resolution  = c(grid_res_x, grid_res_y),
  vals        = NA
)


### ADD VALUES (CONTINENTS)
pos <- cellFromXY(ras, cells_infos[ , c("x", "y")])
ras[][pos] <- 1


footprint2 <- raster(
  xmn         = map_extent["x", "min"],
  xmx         = map_extent["x", "max"],
  ymn         = map_extent["y", "min"],
  ymx         = map_extent["y", "max"],
  crs         = proj4string(map_shp),
  resolution  = c(grid_res_x, grid_res_y),
  vals        = footprint@values
)

crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

r <- raster()
bb <- extent(map_extent["x", "min"], map_extent["x", "max"], map_extent["y", "min"], map_extent["y", "max"])
extent(footprint) <- bb
footprint2 <- setExtent(footprint, bb, keepres=TRUE)

new.rasterA = projectRaster(footprint2, ras) #define the projection and extent

r.stack = stack(new.rasterA, ras)

new_2013 <- crop(extend(new.rasterA, ras), ras)
all.equal(extent(raster_2015), extent(new_2013))