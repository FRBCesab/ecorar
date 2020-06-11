#' Panel of Six World Maps
#'
#' This script produces the Loiseau, Mouquet et al.'s 2020 paper Figure 4, i.e.
#' a panel of six world maps (mammals on the left and birds on the right)
#' representing species richness.
#'
#' @author Nicolas Casajus, \email{nicolas.casajus@@fondationbiodiversite.fr},
#'         Nicolas Loiseau, \email{nicolas.loiseau1@@gmail.com},
#'         Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr}
#'
#' @date 2020/06/11


owarn <- options()$warn
options(warn = -1)


## Parameters ----

let <- 1


## Init Export Device ----

png(
  file      = here::here("figures", paste0(figname, ".png")),
  width     = 24.00,
  height    = 15.75,
  units     = "in",
  res       = 600,
  pointsize = 18
)


## Define Graphical Parameters ----

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
  mfcol    = c(3, 2)
)



for (taxa in taxas) {


  ## Import Data to Map ----

  rasters   <- readRDS(here::here("data", paste0(taxa, "_funk.rds")))

  
  ## Remove Capsian Sea Cells ----
  
  caspian_cells <- c(
    28245, 28246, 28247, 28248, 28938, 28939, 28940, 28941, 29632, 29633, 30326, 
    30327, 30328, 31021, 31022, 31023, 31716, 31717, 31718, 31719, 32411, 32412, 
    32413, 32414, 33107, 33108, 33109, 33110, 33803, 33804, 33805, 33806, 34499, 
    34500, 34501, 34502, 35196, 35197, 35198, 35199, 35892, 35893, 35894, 36585, 
    36586, 36587, 36588, 36589, 37279, 37280, 37281, 37282, 37283, 37284, 37974, 
    37975, 37976, 37977, 37978, 37979, 37981, 38669, 38670, 38671, 38672, 38673, 
    38674, 38675, 38676, 39364, 39365, 39366, 39367, 39368, 39369, 39370, 39371, 
    40061, 40062, 40063, 40064, 40065, 40066
  )

  rasters[][caspian_cells, ] <- NA



  for (j in 1:length(vars_richness)) {



    ## Plot Ocean Background ----

    plot(border, border = NA, col = color_ocean)


    ## Add Graticules ----

    grd <- addGraticules(add = TRUE)


    ## Add Axes ----

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

    
    ## Add Data ----

    rvb <- plotRVB(
      x       = subset(rasters, paste(taxa, vars_richness[j], sep = "_")),
      reverse = FALSE,
      add     = TRUE
    )


    ## Add Figure Box ----

    par(xpd = TRUE)

    plot(border, border = par()$col, lwd = 4, add = TRUE)
    plot(border, border = "white", lwd = 2, add = TRUE)
    
    par(xpd = FALSE)


    ## Add Inset ----

    if (vars_richness[j] == "D75R75") {

      if (taxa == "mammals") {

        
        ## Central America ----

        addInset(rvb, region = c(xmn = -8879810, xmx = -5516986, 
                                 ymn =   632369,  ymx = 3118864),
                 where = c(xmin = -15750000, ymin = -4000000),
                 zoom = 2, title = "", yat = 0.1, case = 1)

        
        ## Madagascar ----

        addInset(rvb, region = c(xmn =  3980550, xmx = 5215893,
                                 ymn = -3517876, ymx = -997417),
                 where = c(xmin = -3000000, ymin = -7000000),
                 zoom = 2.65, title = "", yat = 0.075, case = 3)

        
        ## Asia ----

        addInset(rvb, region = c(xmn =  9047080, xmx = 15780354, 
                                 ymn = -1573001, ymx =   903988),
                 where = c(xmin = 5900000, ymin = -7000000),
                 zoom = 1.65, title = "", yat = 0.12, case = 2)
      }


      if (taxa == "birds") {

       
        ## South America ----
        
        addInset(rvb, region = c(xmn = -8500000, xmx = -6500000, 
                                 ymn = -1300000, ymx =  2000000),
                 where = c(xmin = -14000000, ymin = -4000000),
                 zoom = 2, title = "", yat = 0.1, case = 3)
        
        
        ## Jah'maica (Rastafari) ----
        
        addInset(rvb, region = c(xmn = -7700000, xmx = -5750000, 
                                 ymn =  1600000, ymx =  2750000),
                 where = c(xmin = -5600000, ymin = 2700000),
                 zoom = 2, title = "", yat = 0.1, case = 4)


        ## Asia ----
        
        addInset(rvb, region = c(xmn =  9047080, xmx = 15780354, 
                                 ymn = -1573001, ymx =   903988),
                 where = c(xmin = 5900000, ymin = -7000000),
                 zoom = 1.65, title = "", yat = 0.12, case = 2)
      }
    }


    ## Add Legend Colors ----

    ybottom <-  -6500000
    xstart  <- -16750000
    xleft   <- xstart

    for (i in 1:length(color_richness)) {

      rect(xleft, ybottom, xleft + 35000, ybottom + 350000, border = NA, 
           col = color_richness[i])
      xleft   <- xleft + 35000
    }


    ## Add Legend Box ----

    rect(xstart, ybottom, xleft, ybottom + 350000, border = "white", col = NA)


    ## Add Legend Labels ----

    x <- c(xstart, xleft, xstart + ((xleft - xstart) / 2))
    y <- rep(ybottom, length(x))
    labels <- c(
      0, 
      round(
        max(
          subset(rasters, paste(taxa, vars_richness[j], sep = "_"))[], na.rm = TRUE
        )
      ),
      round(
        max(
          subset(rasters, paste(taxa, vars_richness[j], sep = "_"))[], na.rm = TRUE
        ) / 2
      )
    )

    for (i in 1:length(x)) {
      
      text(x[i], y[i], labels[i], pos = 1, col = par()$col.axis, font = 2, cex = 1.45)  
    }
    
    
    ## Add Legend Title ----

    if (vars_richness[j] ==  "TD_sp") {
      
      text(
        x       = xstart + ((xleft - xstart) / 2),
        y       = ybottom + 350000,
        labels  = paste0("Total number of\n", tolower(taxa), " species"),
        pos     = 3,
        col     = par()$col.axis,
        font    = 2,
        cex     = 1.65
      )
    
    } else {
      
      text(
        x       = xstart + ((xleft - xstart) / 2),
        y       = ybottom + 350000,
        labels  = paste0(
          "Number of ",
          ifelse(vars_richness[j] == "D75R75",
                 "ecologically\nrare ",
                 "ecologically\ncommon "),
          tolower(taxa),
          " species"
        ),
        pos     = 3,
        col     = par()$col.axis,
        font    = 2,
        cex     = 1.65
      )
    }
    

    ## Add Sub-plot ID ----

    text(-16500000, 4900000, letters[let], col = par()$col, font = 2, cex = 2.5)


    ## Add Silhouette ----

    if (vars_richness[j] == "TD_sp") {
      addPhylopic(
        img    = icons[[taxa]],
        x      = .10,
        y      = .82,
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

options(warn = owarn)


## Message ----

usethis::ui_done(
  paste(
    usethis::ui_value(figname), 
    "successfully exported in figures/"
  )
)
