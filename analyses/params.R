#' Set Project Parameters
#'
#' This R script sets project parameters, e.g. main variables, colors values,
#' etc.
#'
#' @author Nicolas Casajus, \email{nicolas.casajus@@fondationbiodiversite.fr},
#'         Nicolas Loiseau, \email{nicolas.loiseau1@@gmail.com},
#'         Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr}
#'
#' @date 2020/06/10


## Define Main Variables ----

taxas         <- c("mammals", "birds")
classes       <- c("D25R25", "AVG", "D75R75")
classes_labs  <- c("Common", "Average", "Rare")
vars_richness <- c("TD_sp", "D75R75", "D25R25")
cc_horizons   <- c("2041_2060", "2061_2080")
iucn_status   <- c("NE", "LC", "TH")
threats_vars  <- c("humanfootprint", "hdi", "conflict", "climate_change",
                   "targetmet_percentagecover")
threats_labs  <- c("Human footprint", "HDI", "Number of conflicts",
                   "Climate change (%)", "Target achievement (%)")


## Define Color Values ----

color_rare    <- "#FF4500"             # ~ Red
color_avg     <- "#00AFBB"             # ~ Turquoise
color_common  <- "#E7B800"             # ~ Orange

color_classes <- c(color_common, color_avg, color_rare)
names(color_classes) <- classes

color_silh    <- "#777777"
light_grey    <- "#888888"
dark_grey     <- "#333333"
par_fg        <- "#666666"
color_ocean   <- "#95D8EB"

color_distinctiveness <- RColorBrewer::brewer.pal(name = "YlGnBu", n = 9)
color_distinctiveness <- colorRampPalette(color_distinctiveness)(255)

color_richness <- RColorBrewer::brewer.pal(name = "YlOrRd", n = 9)
color_richness <- colorRampPalette(color_richness)(255)
color_richness <- c("#AAAAAA", color_richness)

color_iucn_bg <- c("#AAAAAA", "#026666", "#C53131")
color_iucn_fg <- c("#333333", "#F7F7F7", "#F6C8C8")
names(color_iucn_bg) <- names(color_iucn_fg) <- iucn_status

alpha      <- "88"
par_family <- "serif"


gg_theme <- theme_bw() +
  theme(text               = element_text(par_family, colour = dark_grey), 
        axis.text          = element_text(size = 24, colour = dark_grey),
        axis.title         = element_text(size = 24, face = 2),
        legend.title       = element_text(size = 24, face = 2),
        legend.text        = element_text(size = 24),
        legend.background  = element_rect(fill = "white", color = "white"),
        legend.margin      = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

## Define Figure 1 Parameters ----

#...


## Define Figure 2 Parameters ----

pcoa_axes  <- list(mammals = c(2, 3), birds = c(2, 4))
jitter_val <- 500


## Define Figure 3 Parameters ----

#...


## Define Figure 4 Parameters ----

proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

border <- as(raster::extent(c(-180, 180, -90, 90)), "SpatialPolygons")
sp::proj4string(border) <- proj4


## Define Figure 5 Parameters ----

jitter_fac <- c(15, 7, 5)
names(jitter_fac) <- classes
