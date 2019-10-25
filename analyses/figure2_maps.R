#' --------------------------------------------------------------------------   @Header
#'
#' @title Figure 1 - Panel of Pcoa scatterplots
#'
#' @description
#' This R script...
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @date 2019/08/05
#'
#' --------------------------------------------------------------------------   @VVVVVV



#'  -------------------------------------------------------------------------   @Parameters


let <- 1



#'  -------------------------------------------------------------------------   @InitExportDevice


png(
  file      = file.path(path_figs, paste0(figname2, ".png")),
  width     = 12.00 * 2,
  height    =  5.25 * 3,
  units     = "in",
  res       = 600,
  pointsize = 18
)



#'  -------------------------------------------------------------------------   @DefineGraphicalParameters


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



#'  -------------------------------------------------------------------------   @ImportData2Map


  rasters   <- readRDS(file.path(path_data, paste0(taxa, "_funk.rds")))



  for (j in 1:length(vars_richness)) {



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


    rvb <- plotRVB(
      x       = subset(rasters, paste(taxa, vars_richness[j], sep = "_")),
      reverse = FALSE,
      add     = TRUE
    )



#'  -------------------------------------------------------------------------   @AddFigureBox


    par(xpd = TRUE)
    plot(border, border = par()$col, lwd = 4, add = TRUE)
    plot(border, border = "white", lwd = 2, add = TRUE)
    par(xpd = FALSE)



#'  -------------------------------------------------------------------------   @AddInset


    if (vars_richness[j] == "D75R75") {



      if (taxa == "mammals") {



#'  -------------------------------------------------------------------------   @CentralAmerica


        addInset(
          x       = rvb,
          region  = c(
            xmn = -8879810, xmx = -5516986, ymn = 632369, ymx = 3118864
          ),
          where   = c(
            xmin = -15750000, ymin = -5250000
          ),
          zoom    = 2,
          title   = "",
          yat     = 0.1,
          case    = 1
        )



#'  -------------------------------------------------------------------------   @Madagascar


        addInset(
          x       = rvb,
          region  = c(
            xmn = 3980550, xmx = 5215893, ymn = -3517876, ymx = -997417
          ),
          where   = c(
            xmin = -3000000, ymin = -7000000
          ),
          zoom    = 2.65,
          title   = "",
          yat     = 0.075,
          case    = 3
        )




#'  -------------------------------------------------------------------------   @Asia


        addInset(
          x       = rvb,
          region  = c(
            xmn = 9047080, xmx = 15780354, ymn = -1573001, ymx = 903988
          ),
          where   = c(
            xmin = 5900000, ymin = -7000000
          ),
          zoom    = 1.65,
          title   = "",
          yat     = 0.12,
          case    = 2
        )

      }



      if (taxa == "birds") {




#'  -------------------------------------------------------------------------   @SouthAmerica


        addInset(
          x       = rvb,
          region  = c(
            xmn = -8500000, xmx = -6500000, ymn = -1300000, ymx = 2000000
          ),
          where   = c(
            xmin = -14000000, ymin = -5250000
          ),
          zoom    = 2,
          title   = "",
          yat     = 0.1,
          case    = 3
        )



#'  -------------------------------------------------------------------------   @Jamaica


        addInset(
          x       = rvb,
          region  = c(
            xmn = -7700000, xmx = -5750000, ymn = 1600000, ymx = 2750000
          ),
          where   = c(
            xmin = -5600000, ymin = 2700000
          ),
          zoom    = 2,
          title   = "",
          yat     = 0.1,
          case    = 4
        )



#'  -------------------------------------------------------------------------   @Asia


        addInset(
          x       = rvb,
          region  = c(
            xmn = 9047080, xmx = 15780354, ymn = -1573001, ymx = 903988
          ),
          where   = c(
            xmin = 5900000, ymin = -7000000
          ),
          zoom    = 1.65,
          title   = "",
          yat     = 0.12,
          case    = 2
        )
      }
    }



#'  -------------------------------------------------------------------------   @AddColorsLegend


    ybottom <- -6550000
    xstart  <- -16500000
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
      labels  = round(
        max(
          subset(rasters, paste(taxa, vars_richness[j], sep = "_"))[], na.rm = TRUE
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
          subset(rasters, paste(taxa, vars_richness[j], sep = "_"))[], na.rm = TRUE
        ) / 2
      ),
      pos     = 1,
      col     = par()$col.axis,
      font    = 2
    )



#'  -------------------------------------------------------------------------   @AddLegendTitle


    text(
      x       = xstart + ((xleft - xstart) / 2),
      y       = ybottom + 350000,
      labels  = paste0(
        "NUMBER OF ",
        ifelse(
          vars_richness[j] == "TD_sp",
          "TOTAL ",
          ifelse(
            vars_richness[j] == "D75R75",
            "RARE ",
            "COMMON "
          )
        ),
        toupper(taxa),
        " SPECIES"
      ),
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


    if (vars_richness[j] == "TD_sp") {
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
