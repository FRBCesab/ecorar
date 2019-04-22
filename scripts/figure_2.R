################################################################################
###                                                                          ###
###                   FIGURE 2 (PANEL OF 6 MAPS w/ INSETS)                   ###
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
### [1] RColorBrewer_1.1-2 raster_2.8-19      rgeos_0.4-2       rgdal_1.4-3  ###
### [5] sp_1.3-1                                                             ###
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
library(RColorBrewer)
library(rphylopic)


### PARAMETERS -----------------------------------------------------------------

taxas <- c("mammals", "birds")
infos <- c("TD_sp", "D75R75", "D25R25")


### PATH TO DATA ---------------------------------------------------------------

path_data <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/data/"           ### !!!
path_R    <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/RALLL/scripts/"  ### !!!


### LOAD ADDING FUNCTIONS ------------------------------------------------------

source(paste0(path_R, "__addGraticules.R"))
source(paste0(path_R, "__compassRose.R"))
source(paste0(path_R, "__addInset.R"))
source(paste0(path_R, "__plotRVB.R"))
source(paste0(path_R, "__recolorPhylopic.R"))
source(paste0(path_R, "__addPhylopic.R"))


### IMPORT BACKGROUND LAYERS ---------------------------------------------------

study      <- raster(paste0(path_data, "reference_grid_50km.tif"))
# antarctica <- ...


### DEFINE MAP EXTENT (WITH ANTARCTICA) ----------------------------------------

border              <- as(extent(c(-180, 180, -90, 90)), "SpatialPolygons")
proj4string(border) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
border              <- spTransform(border, proj4string(study))


### IMPORT ANIMALS ICONS -------------------------------------------------------

icons <- list()
icons[["mammals"]] <- image_data("5a5dafa2-6388-43b8-a15a-4fd21cd17594", size = 512)[[1]]
icons[["mammals"]] <- icons[["mammals"]][ , dim(icons[["mammals"]])[2]:1, ]
icons[["birds"]]   <- image_data("c3c19b65-cb8d-40bd-b1a6-82a3277bcd4f", size = 512)[[1]]
# icons <- readRDS(paste0(path_data, "icons.rds"))


### INIT EXPORT DEVICE ---------------------------------------------------------

png(
  file       = paste0("~/Desktop/figure-2.png"),
  width      = 12.00 * 2,
  height     =  5.25 * 3,
  units      = "in",
  res        = 600,
  pointsize  = 18
)


### DEFINE GRAPHICAL PARAMETERS ------------------------------------------------

par(
  xaxs      = "i",
  yaxs      = "i",
  family    = "serif",
  mar       = rep(1, 4),
  cex.axis  = 1.25,
  mgp       = c(2, .5, 0),
  tcl       = -0.25,
  xpd       = FALSE,
  new       = FALSE,
  fig       = c(0, 1, 0, 1),
  col       = "#666666",
  col.axis  = "#666666",
  fg        = "#666666",
  mfcol     = c(3, 2)
)


let <- 1

for (k in 1:length(taxas)) {


  ### IMPORT DATA --------------------------------------------------------------

  datas   <- readRDS(paste0(path_data, "funk_", taxas[k], ".rds"))


  for (j in 1:length(infos)) {


    ### DEFINE COLORS GRADIENT -------------------------------------------------

    cols <- c("#aaaaaa", colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(255))


    ### DEFINE MAP EXTENT ------------------------------------------------------

    plot(border, border = NA, col = "#95D8EB")


    ### ADD GRATICULES AND AXIS ------------------------------------------------

    grd <- addGraticules(add = TRUE, line.color = "#aaaaaa")

    par(mgp = c(.5, -.1, 0))
    axis(
      side      = 1,
      at        = grd[[1]][ , "x"],
      labels    = paste0(grd[[1]][ , "label"], "°", grd[[1]][ , "direction"]),
      cex.axis  = .75,
      lwd       = 0
    )

    par(mgp = c(.5, -.45, 0))
    axis(
      side      = 4,
      at        = grd[[2]][ , "y"],
      labels    = paste0(grd[[2]][ , "label"], "°", grd[[2]][ , "direction"]),
      cex.axis  = .75,
      lwd       = 0
    )

    par(mgp = c(.5, 0, 0))
    axis(
      side      = 3,
      at        = grd[[1]][ , "x"],
      labels    = paste0(grd[[1]][ , "label"], "°", grd[[1]][ , "direction"]),
      cex.axis  = .75,
      lwd       = 0
    )

    par(mgp = c(.5, -.35, 0))
    axis(
      side      = 2,
      at        = grd[[2]][ , "y"],
      labels    = paste0(grd[[2]][ , "label"], "°", grd[[2]][ , "direction"]),
      cex.axis  = .75,
      lwd       = 0
    )


    ### ADD RASTER DATA --------------------------------------------------------

    rvb <- plotRVB(
      x        = subset(datas, paste(taxas[k], infos[j], sep = "_")),
      reverse  = FALSE,
      add      = TRUE
    )


    ### ADD FIGURE BOX  --------------------------------------------------------

    par(xpd = TRUE)
    plot(border, border = par()$col, lwd = 4, add = TRUE)
    plot(border, border = "white", lwd = 2, add = TRUE)
    par(xpd = FALSE)


    ### ADD MAP INSETS  --------------------------------------------------------

    if (infos[j] == "D75R75") {

      if (taxas[k] == "mammals") {

        addInset(
          x       = rvb,
          region  = c(xmn = -8879810, xmx = -5516986, ymn = 632369, ymx = 3118864),
          where   = c(xmin = -15750000, ymin = -5250000),
          zoom    = 2,
          title   = "",
          yat     = 0.1,
          case    = 1
        )

        addInset(
          x       = rvb,
          region  = c(xmn = 3980550, xmx = 5215893, ymn = -3517876, ymx = -997417),
          where   = c(xmin = -3000000, ymin = -7000000),
          zoom    = 2.65,
          title   = "",
          yat     = 0.075,
          case    = 3
        )

        addInset(
          x       = rvb,
          region  = c(xmn = 9047080, xmx = 15780354, ymn = -1573001, ymx = 903988),
          where   = c(xmin = 5900000, ymin = -7000000),
          zoom    = 1.65,
          title   = "",
          yat     = 0.12,
          case    = 2
        )

        # addInset(
        #   x       = rvb,
        #   region  = c(xmn = 11053993, xmx = 12457308, ymn = 420761, ymx = 3458091),
        #   where   = c(xmin = 14460000, ymin = 1450000),
        #   zoom    = 1.85,
        #   title   = "",
        #   yat     = 0.18,
        #   case    = 4
        # )
      }

      if (taxas[k] == "birds") {

        addInset(
          x       = rvb,
          region  = c(xmn = -8500000, xmx = -6500000, ymn = -1300000, ymx = 2000000),
          where   = c(xmin = -14000000, ymin = -5250000),
          zoom    = 2,
          title   = "",
          yat     = 0.1,
          case    = 3
        )

        addInset(
          x       = rvb,
          region  = c(xmn = -7700000, xmx = -5750000, ymn = 1600000, ymx = 2750000),
          where   = c(xmin = -5600000, ymin = 2700000),
          zoom    = 2,
          title   = "",
          yat     = 0.1,
          case    = 4
        )

        addInset(
          x       = rvb,
          region  = c(xmn = 9047080, xmx = 15780354, ymn = -1573001, ymx = 903988),
          where   = c(xmin = 5900000, ymin = -7000000),
          zoom    = 1.65,
          title   = "",
          yat     = 0.12,
          case    = 2
        )
      }
    }


    ### ADD LEGEND COLORS  -----------------------------------------------------

    xstart  <- -17050000
    xleft   <- xstart
    ybottom <- -6550000

    for (i in 1:length(cols)) {

      rect(
        xleft    = xleft,
        ybottom  = ybottom,
        xright   = xleft + 35000,
        ytop     = ybottom + 350000,
        border   = NA,
        col      = cols[i]
      )
      xleft   <- xleft + 35000
    }

    rect(
      xleft    = xstart,
      ybottom  = ybottom,
      xright   = xleft,
      ytop     = ybottom + 350000,
      border   = "white",
      col      = NA
    )


    ### ADD LEGEND LABELS  -----------------------------------------------------

    text(
      x       = xstart,
      y       = ybottom,
      labels  = 0,
      pos     = 1,
      col     = par()$col.axis,
      font    = 2
    )

    text(
      x       = xleft,
      y       = ybottom,
      labels  = round(
        max(
          subset(datas, paste(taxas[k], infos[j], sep = "_"))[], na.rm = TRUE
        )
      ),
      pos     = 1,
      col     = par()$col.axis,
      font    = 2
    )

    text(
      x       = xstart + ((xleft - xstart) / 2),
      y       = ybottom,
      labels  = round(
        max(
          subset(datas, paste(taxas[k], infos[j], sep = "_"))[], na.rm = TRUE
        ) / 2
      ),
      pos     = 1,
      col     = par()$col.axis,
      font    = 2
    )


    ### ADD LEGEND TITLE  ------------------------------------------------------

    text(
      x       = xstart + ((xleft - xstart) / 2),
      y       = ybottom + 350000,
      labels  = paste0(
        "Number of ",
        taxas[k],
        " species (",
        toupper(gsub("_", " ", infos[j])),
        ")"
      ),
      pos     = 3,
      col     = par()$col.axis,
      font    = 2
    )


    ### ADD NORTH ARROW  -------------------------------------------------------

    compassRose(
      x        = -14000000,
      y        = 4900000,
      cex.cr   = .65,
      cex.let  = .65
    )


    ### ADD SUBPLOT LABEL  -----------------------------------------------------

    text(
      x       = -16500000,
      y       = 4900000,
      labels  = LETTERS[let],
      col     = par()$col,
      font    = 2,
      cex     = 2
    )

    let <- let + 1


    ### ADD ICON ---------------------------------------------------------------

    if (infos[j] == "TD_sp") {
      addPhylopic(
        img    = icons[[taxas[k]]],
        x      = .05,
        y      = .60,
        ysize  = .075,
        alpha  = 1,
        color  = "#777777",
        AR     = ifelse(taxas[k] == "mammals", 1.5, 2)
      )
    }
  }
}


dev.off()
