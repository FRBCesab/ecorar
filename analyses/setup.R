#' --------------------------------------------------------------------------   @Header
#'
#' @title Project setup
#'
#' @description
#' This R script...
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @date 2019/10/25
#'
#' --------------------------------------------------------------------------   @Header



#'  -------------------------------------------------------------------------   @SetOptionsParameters


options(warn = -1)



#' ----------------------------------------------------------------------------- @InstallCranLibs


pkgs <- c(
  "png",
  "ggplot2",
  "cowplot",
  "grid",
  "gridExtra",
  "ape",
  "treeio",
  "tidytree",
  "sp",
  "rgdal",
  "rgeos",
  "raster",
  "RColorBrewer",
  "leaflet",
  "BiocManager"
)

nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)



#' ----------------------------------------------------------------------------- @InstallBioManagerLibs


if (!("ggtree" %in% installed.packages())) { BiocManager::install("ggtree") }



#' ----------------------------------------------------------------------------- @LoadLibs


pkgs <- c(pkgs, "ggtree")
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

if (sum(ip) != length(pkgs)) { cat("Some packages failed to load.\n") }

rm(list = c("pkgs", "nip", "ip"))



#' ----------------------------------------------------------------------------- @LoadRFunctions


rfun <- list.files(path = "R", pattern = "^__.+\\.R$", full.names = TRUE)
rfun <- unlist(lapply(rfun, source, verbose = FALSE))

rm(list = "rfun")



#' ----------------------------------------------------------------------------- @CreateFolders


dir_names <- c(
  "data",
  "figures"
)

dir_vars  <- c(
  "path_data",
  "path_figs"
)

dirs <- lapply(

  X   = 1:length(dir_names),

  FUN = function(i) {

    dir.create(
      path         = dir_names[i],
      showWarnings = FALSE,
      recursive    = TRUE
    )

    assign(
      x     = dir_vars[i],
      value = dir_names[i],
      envir = .GlobalEnv
    )
  }
)

rm(list = c("dir_names", "dir_vars", "dirs"))



#' ----------------------------------------------------------------------------- @SetFigureNames


fig_names <- c(
  "Figure_1",
  "Figure_2",
  "Figure_3",
  "Figure_4"
)

fig_vars  <- c(
  "figname1",
  "figname2",
  "figname3",
  "figname4"
)

figs <- lapply(

  X   = 1:length(fig_names),

  FUN = function(i) {

    assign(
      x     = fig_vars[i],
      value = fig_names[i],
      envir = .GlobalEnv
    )
  }
)

rm(list = c("fig_names", "fig_vars", "figs"))



#'  -------------------------------------------------------------------------   @GlobalParameters


taxas          <- c("mammals", "birds")

classes        <- c("D25R25", "AVG", "D75R75")
classes_labs   <- c("Common", "Average", "Rare")
vars_richness  <- c("TD_sp", "D75R75", "D25R25")

cc_horizons    <- c("2041_2060", "2061_2080")

threats_vars   <- c(
  "humanfootprint", "hdi", "conflict", "climate_change", "targetmet_percentagecover"
)
threats_labs   <- c(
  "Human footprint", "HDI", "Number of conflicts", "Climate change (%)", "Target achievement (%)"
)

iucn_status    <- c("NE", "LC", "TH")



#'  -------------------------------------------------------------------------   @ColorsParameters


color_rare    <- "#ff4500"             # ~ Red
color_avg     <- "#00afbb"             # ~ Turquoise
color_common  <- "#e7b800"             # ~ Orange

color_classes <- c(color_common, color_avg, color_rare)
names(color_classes) <- classes

color_silh    <- "#777777"
light_grey    <- "#888888"
dark_grey     <- "#333333"
par_fg        <- "#666666"
color_ocean   <- "#95d8eb"

color_distinctiveness <- RColorBrewer::brewer.pal(name = "YlGnBu", n = 9)
color_distinctiveness <- colorRampPalette(color_distinctiveness)(255)

color_richness <- RColorBrewer::brewer.pal(name = "YlOrRd", n = 9)
color_richness <- colorRampPalette(color_richness)(255)
color_richness <- c("#aaaaaa", color_richness)

color_iucn_bg <- c("#aaaaaa", "#026666", "#c53131")
color_iucn_fg <- c("#333333", "#f7f7f7", "#f6c8c8")
names(color_iucn_bg) <- names(color_iucn_fg) <- iucn_status

alpha      <- "88"
par_family <- "serif"


#'  -------------------------------------------------------------------------   @Fig1Parameters


pcoa_axes <- list(
  mammals = c(2, 3),
  birds   = c(2, 4)
)

jitter_val <- 500



#'  -------------------------------------------------------------------------   @Fig2Parameters


proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"



#'  -------------------------------------------------------------------------   @Fig4Parameters


jitter_fac <- c(15, 7, 5)
names(jitter_fac) <- classes




#'  -------------------------------------------------------------------------   @LoadSpeciesSilh


icons <- lapply(taxas, function(x) {
    png::readPNG(source = file.path(path_data, paste0(x, "_silhouette.png")))
  }
)

names(icons) <- taxas



#'  -------------------------------------------------------------------------   @ImportGridRaster


study <- raster(file.path(path_data, "reference_grid_50km.tif"))



#'  -------------------------------------------------------------------------   @DefineMapExtent


border              <- as(extent(c(-180, 180, -90, 90)), "SpatialPolygons")
proj4string(border) <- proj4
border              <- spTransform(border, proj4string(study))



#'  -------------------------------------------------------------------------   @LoadData4Fig2


estimated_d <- get(
  load(
    file = file.path(
      path_data,
      "estimated_d_for_figure2.RData"
    )
  )
)



#'  -------------------------------------------------------------------------   @LoadAllData


datas   <- get(load(file = file.path(path_data, "alldata_for_figures.RData")))



#'  -------------------------------------------------------------------------   @PrepareIUCNData

iucn <- datas[!is.na(datas[ , "iucn_status"]), ]
iucn <- iucn[!is.na(iucn[ , "dr_class"]), c("iucn_status", "dr_class", "class")]

iucn <- tapply(
  X     = datas[ , "iucn_status"],
  INDEX = list(
    datas[ , "iucn_status"], datas[ , "dr_class"], datas[ , "class"]
  ),
  FUN   = function(x) { length(x) }
)

iucn <- list(
  mammals = iucn[ , , "mammals"],
  birds   = iucn[ , , "birds"]
)

iucn <- lapply(
  X    = iucn,
  FUN  = function(x) {
    x[is.na(x)] <- 0
    return(x)
  }
)

iucn <- lapply(
  X    = iucn,
  FUN  = function(y) { apply(y, 2, function(x) x / sum(x) * 100) }
)
