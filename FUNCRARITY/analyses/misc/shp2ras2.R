################################################################################
###                                                                          ###
###                            SHAPEFILE TO RASTER                           ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### AUTHORS : The Three Nico                                                 ###
### DATE    : 2019/04/08                                                     ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### > sessionInfo()                                                          ###
###                                                                          ###
### R version 3.5.3 (2019-03-11)                                             ###
### Platform: x86_64-apple-darwin18.2.0 (64-bit)                             ###
### Running under: macOS Mojave 10.14.4                                      ###
###                                                                          ###
### locale:                                                                  ###
### [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8        ###
###                                                                          ###
### attached base packages:                                                  ###
### [1] stats     graphics  grDevices utils     datasets  methods   base     ###
###                                                                          ###
### other attached packages:                                                 ###
### [1] rgdal_1.4-3   rgeos_0.4-2   raster_2.8-19    sp_1.3-1                ###
###                                                                          ###
### loaded via a namespace (and not attached):                               ###
### [1] compiler_3.5.3   Rcpp_1.0.1       codetools_0.2-16    grid_3.5.3     ###
### [5] lattice_0.20-38                                                      ###
###                                                                          ###
################################################################################



rm(list = ls())



library(sp)
library(rgdal)
library(rgeos)
library(raster)



### PATH TO DATA
path_data <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/data/"       ### !!!


### IMPORT STUDY AREA SHAPEFILE (WORLD GRID)
map_shp <- readOGR(
  dsn    = paste0(path_data, "ReferenceGrid50km"),                              ### !!!
  layer  = "gridLand50km"
)


### GET MAP INFOS
map_extent      <- map_shp@bbox
poly_coords     <- lapply(map_shp@polygons, function(x) x@Polygons[[1]]@coords)
grid_res_x      <- unique(unlist(lapply(poly_coords, function(x) max(x[ , 1]) - min(x[ , 1]))))
grid_res_y      <- unique(unlist(lapply(poly_coords, function(x) max(x[ , 2]) - min(x[ , 2]))))

cells_infos <- data.frame(
  cell_id = map_shp@data[ , 1],
  gCentroid(map_shp, byid = TRUE)
)


### CREATE TEMPLATE RASTER (GEOGRAPHICAL INFOS)
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


### EXPORT STUDY AREA RASTER
writeRaster(
  x          = ras,
  filename   = paste0(path_data, "reference_grid_50km.tif"),                     ### !!!
  format     = "GTiff",
  overwrite  = TRUE
)


### CONVERT MAMMALS / BIRDS TO STACKED RASTER

taxa <- c("birds", "mammals")

for (k in 1:length(taxa)) {


  # Import mammals and birds data
  load(paste0(path_data, taxa[k], "/SES_total_", taxa[k], ".RData"))
  datas <- eval(parse(text = paste0("SES_total_", taxa[k])))


  # Get good cells id order
  datas <- merge(cells_infos, datas, by.x = "cell_id", by.y = "cell")
  pos   <- cellFromXY(ras, datas[ , c("x", "y")])


  # Add data to raster (stacked raster)
  sar <- ras
  for (i in 4:ncol(datas)) {
    tmp <- ras ; tmp[] <- NA
    tmp[][pos] <- datas[ , i]
    if (i == 4) {
      sar <- tmp
    } else {
      sar <- stack(sar, tmp)
    }
  }

  # Add layers names
  names(sar) <- paste(taxa[k], colnames(datas)[-c(1:3)], sep = "_")

  # Export stacked layer (RDS file)
  saveRDS(
    object  = sar,
    file    = paste0(path_data, "SES_total_", taxa[k], ".rds")
  )
}
