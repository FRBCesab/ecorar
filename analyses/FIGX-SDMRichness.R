################################################################################
###                                                                          ###
###                   FIGURE 3 (PANEL OF 6 MAPS w/ INSETS)                   ###
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
### R version 3.5.3 (2019-03-11) -- "Great Truth"                            ###
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
### [5] sp_1.3-1   png_0.1-7                                                 ###
###                                                                          ###
################################################################################



rm(list = ls())



#'  -------------------------------------------------------------------------   @Parameters


root      <- "/Users/nicolascasajus/OneDrive/OneDrive - Fondation Biodiversité/MySpace/GROUPS/FREE/01-Loiseau/RALLL/FUNCRARITY/"
source(file.path(root, "graphsParameters.R"))

filename <- "Figure_X"

let <- 1



#'  -------------------------------------------------------------------------   @LoadAddings


library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(RColorBrewer)
library(png)

addings <- list.files(path = path_R, pattern = "\\.R$", full.names = TRUE)
for (i in 1:length(addings)) { source(addings[i]) }



#'  -------------------------------------------------------------------------   @ImportGrid


study <- raster(file.path(path_data, "reference_grid_50km.tif"))



#'  -------------------------------------------------------------------------   @DefineMapExtent


border              <- as(extent(c(-180, 180, -90, 90)), "SpatialPolygons")
proj4string(border) <- proj4
border              <- spTransform(border, proj4string(study))



#'  -------------------------------------------------------------------------   @ImportSpeciesSilhouettes


icons <- list()

for (taxa in taxas) {

  icons[[taxa]] <- readPNG(
      source = file.path(path_data, taxa, paste0(taxa, "_silhouette.png"))
  )
}



#'  -------------------------------------------------------------------------   @InitExportDevice


png(
  file       = file.path(path_figs, paste0(filename, "_Richness.png")),
  width      = 12.00 * 2,
  height     =  5.25 * 4,
  units      = "in",
  res        = 600,
  pointsize  = 18
)



#'  -------------------------------------------------------------------------   @DefineGraphicalParameters


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
  mfcol     = c(4, 2)
)


varnames <- c("observations", "1979-2013", "2041-2060", "2061-2080")
zmaxs <- c(420, 1800)
names(zmaxs) <- taxas


for (k in 1:length(taxas)) {



  taxa <- taxas[k]



  for (j in 1:length(varnames)) {


    varname <- varnames[j]



#'  -------------------------------------------------------------------------   @ImportData2Map


    datas <- get(
      load(
        file.path(
          path_data,
          taxa,
          paste0("biomod-richness_", taxa, "_", varname, ".RData")
        )
      )
    )



#'  -------------------------------------------------------------------------   @PlotOcean


    plot(border, border = NA, col = color_ocean)



#'  -------------------------------------------------------------------------   @AddGraticules


    grd <- addGraticules(add = TRUE)



#'  -------------------------------------------------------------------------   @AddAxes


    par(mgp = c(0.5, -0.1, 0))
    axis(
      side      = 1,
      at        = grd[[1]][ , "x"],
      labels    = paste0(grd[[1]][ , "label"], "°", grd[[1]][ , "direction"]),
      cex.axis  = 0.75,
      lwd       = 0
    )

    par(mgp = c(0.5, -0.45, 0))
    axis(
      side      = 4,
      at        = grd[[2]][ , "y"],
      labels    = paste0(grd[[2]][ , "label"], "°", grd[[2]][ , "direction"]),
      cex.axis  = 0.75,
      lwd       = 0
    )

    par(mgp = c(0.5, 0, 0))
    axis(
      side      = 3,
      at        = grd[[1]][ , "x"],
      labels    = paste0(grd[[1]][ , "label"], "°", grd[[1]][ , "direction"]),
      cex.axis  = 0.75,
      lwd       = 0
    )

    par(mgp = c(0.5, -0.35, 0))
    axis(
      side      = 2,
      at        = grd[[2]][ , "y"],
      labels    = paste0(grd[[2]][ , "label"], "°", grd[[2]][ , "direction"]),
      cex.axis  = 0.75,
      lwd       = 0
    )



#'  -------------------------------------------------------------------------   @AddData


    rvb <- plotRVB(
      x        = datas,
      reverse  = FALSE,
      add      = TRUE,
      zmax     = zmaxs[taxa]
    )



#'  -------------------------------------------------------------------------   @AddFigureBox


    par(xpd = TRUE)
    plot(border, border = par()$col, lwd = 4, add = TRUE)
    plot(border, border = "white", lwd = 2, add = TRUE)
    par(xpd = FALSE)



#'  -------------------------------------------------------------------------   @AddColorsLegend


    ybottom <- -6550000
    xstart  <- -16500000 # -17050000
    xleft   <- xstart

    for (i in 1:length(color_richness)) {

      rect(
        xleft    = xleft,
        ybottom  = ybottom,
        xright   = xleft + 35000,
        ytop     = ybottom + 350000,
        border   = NA,
        col      = color_richness[i]
      )

      xleft   <- xleft + 35000
    }



#'  -------------------------------------------------------------------------   @AddLegendBox


    rect(
      xleft    = xstart,
      ybottom  = ybottom,
      xright   = xleft,
      ytop     = ybottom + 350000,
      border   = "white",
      col      = NA
    )



#'  -------------------------------------------------------------------------   @AddLegendLabel


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
      labels  = zmaxs[taxa],
      pos     = 1,
      col     = par()$col.axis,
      font    = 2
    )

    text(
      x       = xstart + ((xleft - xstart) / 2),
      y       = ybottom,
      labels  = zmaxs[taxa] / 2,
      pos     = 1,
      col     = par()$col.axis,
      font    = 2
    )



#'  -------------------------------------------------------------------------   @AddLegendTitle


    text(
      x       = xstart + ((xleft - xstart) / 2),
      y       = ybottom + 350000,
      labels  = paste0("Species richness (", varname, ")"),
      pos     = 3,
      col     = par()$col.axis,
      font    = 2
    )



#'  -------------------------------------------------------------------------   @AddSubplotID


    text(
      x       = -16500000,
      y       = 4900000,
      labels  = letters[let],
      col     = par()$col,
      font    = 2,
      cex     = 2
    )



#'  -------------------------------------------------------------------------   @AddSilhouette


    if (j == 1) {
      addPhylopic(
        img    = icons[[taxa]],
        x      = .05,
        y      = .60,
        ysize  = .075,
        alpha  = 1,
        color  = "#777777",
        AR     = ifelse(taxa == "mammals", 1.5, 2)
      )
    }


    let <- let + 1

  }
}


dev.off()
