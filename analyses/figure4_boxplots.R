#' --------------------------------------------------------------------------   @Header
#'
#' @title Figure 4 - Panel of Threats Boxplots
#'
#' @description
#' This R script...
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr},
#'         Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @date 2019/08/05
#'
#' --------------------------------------------------------------------------   @VVVVVV




for (horizon in cc_horizons) {



#'  -------------------------------------------------------------------------   @InitExportDevice


  png(
    file      = file.path(path_figs, paste0(figname4, "_", horizon, ".png")),
    width     = 12.00,
    height    =  7.00,
    units     = "in",
    res       = 600,
    pointsize = 18
  )



#'  -------------------------------------------------------------------------   @DefineGraphicalParameters


  par(
    family   = par_family,
    mar      = rep(1, 4),
    mgp      = c(2, .5, 0),
    tcl      = -0.25,
    xpd      = FALSE,
    new      = FALSE,
    fig      = c(0, 1, 0, 1),
    fg       = par_fg,
    col      = par_fg,
    col.axis = par_fg,
    cex.axis = 0.75
  )

  par(mfcol = c(2, 5))



  for (taxa in taxas) {



#'  -------------------------------------------------------------------------   @SubplotMargins


    if (taxa == "mammals") {

      par(mar = c(0.5, 1.0, 1.0, 0.5))

    }

    if (taxa == "birds") {

      par(mar = c(1.0, 1.0, 0.5, 0.5))

    }



#'  -------------------------------------------------------------------------   @EmptyPlotIUCN


    plot(0,
      xlim = c(-10, 110),
      ylim = c(0.5, 3.5),
      axes = FALSE,
      ann  = FALSE,
      type = "n"
    )



#'  -------------------------------------------------------------------------   @AddBackground


    grid()
    box()



#'  -------------------------------------------------------------------------   @AddIUCNStickers


    for (classe in classes) {

      for (status in iucn_status) {

        points(
          x   = iucn[[taxa]][status, classe],
          y   = which(classes == classe),
          cex = 4,
          col = color_iucn_bg[status],
          pch = 19
        )

        text(
          x      = iucn[[taxa]][status, classe],
          y      = which(classes == classe),
          labels = status,
          cex    = 0.85,
          col    = color_iucn_fg[status],
          family = "sans",
          font   = 2
        )

      } # e_o iucn_status
    } # e_o classes



#'  -------------------------------------------------------------------------   @AddxAxis


    if (taxa == "birds") {

      par(mgp = c(0.5, -0.05, 0))
      axis(side = 1, lwd = 0, font = 2)

    }



#'  -------------------------------------------------------------------------   @AddyAxis


    par(mgp = c(0.5, 0, 0))
    axis(side = 2, at = 1:3, labels = classes_labs, lwd = 0, font = 2)



#'  -------------------------------------------------------------------------   @AddTitle


    if (taxa == "mammals") {

      mtext(text = "IUCN status (%)", side = 3, font = 2, cex = 0.85)

    }



#'  -------------------------------------------------------------------------   @AddSilhouette


    addPhylopic(
      img   = icons[[taxa]],
      x     = 0.85,
      y     = 0.92,
      ysize = 0.20,
      alpha = 1,
      color = color_silh,
      AR    = 0.5
    )

  } # e_o taxas



#'  -------------------------------------------------------------------------   @ThreatsPlots


  comp <- 3



  for (threat in threats_vars) {



    if (threat == "climate_change") {

      colname <- paste0("cc_", horizon)

    } else {

      colname <- threat

    }

    subdatas <- datas[!is.na(datas[ , colname]), c("class", "dr_class", colname)]



    for (taxa in taxas) {



      if (comp / 2 == round(comp / 2)) {

        par(mar = c(1.0, 1.0, 0.5, 0.5))

      } else {

        par(mar = c(0.5, 1.0, 1.0, 0.5))

      }

      ssubdatas <- subdatas[subdatas[ , "class"] == taxa, ]



#'  -------------------------------------------------------------------------   @EmptyPlot


      xmin <- min(subdatas[ , colname])
      xmax <- max(subdatas[ , colname])

      xmax <- xmax + (xmax - xmin) * 0.15

      plot(0,
        xlim = c(xmin, xmax),
        ylim = c(0.5, 3.5),
        axes = FALSE,
        ann  = FALSE,
        type = "n"
      )



#'  -------------------------------------------------------------------------   @AddBackground


      grid()
      box()



      for (classe in classes){



#'  -------------------------------------------------------------------------   @AddJitterPoints


        x_coords <- ssubdatas[ssubdatas[ , "dr_class"] == classe, colname]
        x_coords <- x_coords[!is.na(x_coords)]

        y_coords <- jitter(
          x      = rep(
            x     = which(classes == classe),
            times = length(x_coords)
          ),
          factor = jitter_fac[classe]
        )

        points(
          x   = x_coords,
          y   = y_coords,
          pch = 21,
          bg  = paste0(color_classes[classe], alpha),
          col = "transparent",
          cex = 0.5
        )



#'  -------------------------------------------------------------------------   @AddBoxplot


        boxplot(
          x          = x_coords,
        	names      = NULL,
        	horizontal = TRUE,
        	col        = "transparent",
          border     = "black",
        	lty        = 1,
        	lwd        = 1,
        	at         = which(classes == classe),
        	add        = TRUE,
          outline    = FALSE,
          axes       = FALSE,
          width      = 1.8
        )



#'  -------------------------------------------------------------------------   @AddVerticalLineCC


        if (comp > 8) {

          abline(v = 0, lwd = 1, lty = 2)

        }



#'  -------------------------------------------------------------------------   @Addpvalues


        pval <- pvalues[
          which(
            pvalues[ , "taxa"]     == taxa &
            pvalues[ , "threat"]   == threat &
            pvalues[ , "dr_class"] == classe
          ), "signif"
        ]

        if (length(pval) > 0) {

          text(
            x      = max(x_coords),
            y      = which(classes == classe),
            labels = pval,
            cex    = 0.85,
            font   = 2,
            pos    = 4
          )

        }

      } # e_o classes



#'  -------------------------------------------------------------------------   @AddxAxis


      if (comp / 2 == round(comp / 2)) {

        par(mgp = c(0.5, -0.05, 0))
        axis(side = 1, lwd = 0, font = 2)

      }



#'  -------------------------------------------------------------------------   @AddTitle


      if (comp / 2 != round(comp / 2)) {

        mtext(
          text = threats_labs[which(threats_vars == threat)],
          side = 3,
          font = 2,
          cex  = 0.85
        )
      }

      comp <- comp + 1

    } # e_o taxas
  } # e_o threats_vars


  dev.off()

} # e_o cc_horizons



#'  -------------------------------------------------------------------------   @ResetOptionsParameters


options(warn = 0)
