#' Import Required Datasets
#'
#' This R script imports required data to produce all figures.
#'
#' @author Nicolas Casajus, \email{nicolas.casajus@@fondationbiodiversite.fr},
#'         Nicolas Loiseau, \email{nicolas.loiseau1@@gmail.com},
#'         Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr}
#'
#' @date 2020/06/10


## Load All Data ----

datas <- get(load(here::here("data", "alldata_for_figures.RData")))


## Load Species Silhouettes ----

icons <- lapply(taxas, function(x)
    png::readPNG(here::here("data", "png", paste0(x, "_silhouette.png")))
)
names(icons) <- taxas


## Import Grid Raster ----

study <- raster::raster(here::here("data", "reference_grid_50km.tif"))


## Project Map Extent ----

border <- sp::spTransform(border, sp::proj4string(study))


## Load Specific Data (Phylogeny Insets) ----

estimated_d <- get(load(here::here("data", "estimated_d_for_figure2.RData")))


## Clean IUCN Data (Violin Plots) ----

iucn <- datas[!is.na(datas$"iucn_status"), ]
iucn <- iucn[!is.na(iucn$"dr_class"), c("iucn_status", "dr_class", "class")]

iucn <- tapply(datas$"iucn_status",
               list(datas$"iucn_status", datas$"dr_class", datas$"class"),
               function(x) length(x))

iucn <- list(mammals = iucn[ , , "mammals"], birds = iucn[ , , "birds"])

iucn <- lapply(iucn, function(x) {
               x[is.na(x)] <- 0
               return(x)})

iucn <- lapply(iucn, function(y) apply(y, 2, function(x) x / sum(x) * 100))
