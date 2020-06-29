png(
  file      = here::here("figures", paste0(figname, ".png")),
  width     = 12.00 * 1,
  height    =  5.25 * 2,
  units     = "in",
  res       = 600,
  pointsize = 18
)

par(
  xaxs     = "i",
  yaxs     = "i",
  family   = par_family,
  mar      = rep(1, 4),
  cex.axis = 1.25,
  mgp      = c(2, .5, 0),
  tcl      = -0.25,
  xpd      = FALSE,
  new      = FALSE,
  fig      = c(0, 1, 0, 1),
  col      = par_fg,
  col.axis = par_fg,
  fg       = par_fg,
  mfcol    = c(2, 1)
)



for (taxa in taxas) {



#'  -------------------------------------------------------------------------   @ImportData2Map


  rasters   <- raster::raster(here::here("data", paste0(taxa, "_Null_D75R75.tif")))


  caspian_cells <- c(
    28245, 28246, 28247, 28248, 28938, 28939, 28940, 28941, 29632, 29633, 30326, 30327,
    30328, 31021, 31022, 31023, 31716, 31717, 31718, 31719, 32411, 32412, 32413, 32414,
    33107, 33108, 33109, 33110, 33803, 33804, 33805, 33806, 34499, 34500, 34501, 34502,
    35196, 35197, 35198, 35199, 35892, 35893, 35894, 36585, 36586, 36587, 36588, 36589,
    37279, 37280, 37281, 37282, 37283, 37284, 37974, 37975, 37976, 37977, 37978, 37979,
    37981, 38669, 38670, 38671, 38672, 38673, 38674, 38675, 38676, 39364, 39365, 39366,
    39367, 39368, 39369, 39370, 39371, 40061, 40062, 40063, 40064, 40065, 40066
  )

  rasters[][caspian_cells] <- NA


#'  -------------------------------------------------------------------------   @PlotOcean


  plot(border, border = NA, col = color_ocean)



#'  -------------------------------------------------------------------------   @AddGraticules


  grd <- addGraticules(add = TRUE)



#'  -------------------------------------------------------------------------   @AddAxes


  par(mgp = c(0.5, -0.1, 0))
  axis(
    side     = 1,
    at       = grd[[1]][ , "x"],
    labels   = paste0(grd[[1]][ , "label"], "°", grd[[1]][ , "direction"]),
    cex.axis = 0.75,
    lwd      = 0
  )

  par(mgp = c(0.5, -0.45, 0))
  axis(
    side     = 4,
    at       = grd[[2]][ , "y"],
    labels   = paste0(grd[[2]][ , "label"], "°", grd[[2]][ , "direction"]),
    cex.axis = 0.75,
    lwd      = 0
  )

  par(mgp = c(0.5, 0, 0))
  axis(
    side     = 3,
    at       = grd[[1]][ , "x"],
    labels   = paste0(grd[[1]][ , "label"], "°", grd[[1]][ , "direction"]),
    cex.axis = 0.75,
    lwd      = 0
  )

  par(mgp = c(0.5, -0.35, 0))
  axis(
    side     = 2,
    at       = grd[[2]][ , "y"],
    labels   = paste0(grd[[2]][ , "label"], "°", grd[[2]][ , "direction"]),
    cex.axis = 0.75,
    lwd      = 0
  )



#'  -------------------------------------------------------------------------   @AddData


  if (taxa == "mammals") breaks   <- c(-1.96, 1.96)
  if (taxa == "birds")   breaks   <- c(0, 0.4712917)
  
  rvb <- plotRVB2(
    x         = rasters,
    n_classes = 255,
    breaks    = breaks,
    palettes  = list(c("#003c30", "#80cdc1"), "#ffffff", c("#dfc27d", "#543005")),
    add       = TRUE
  )



#'  -------------------------------------------------------------------------   @AddFigureBox


  par(xpd = TRUE)
  plot(border, border = par()$col, lwd = 4, add = TRUE)
  plot(border, border = "white", lwd = 2, add = TRUE)
  par(xpd = FALSE)





#'  -------------------------------------------------------------------------   @AddLegendTitle


#'  -------------------------------------------------------------------------   @AddSilhouette


  addPhylopic(
    img    = icons[[taxa]],
    x      = .075,
    y      = .80,
    ysize  = .065,
    alpha  = 1,
    color  = "#777777",
    AR     = ifelse(taxa == "mammals", 1.5, 2)
  )
}

dev.off()
