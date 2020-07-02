#' Panel of Violin Plots
#'
#' This script produces the Loiseau, Mouquet et al.'s 2020 paper Figure 5, i.e.
#' a panel of several violin plots.
#'
#' @author Nicolas Casajus, \email{nicolas.casajus@@fondationbiodiversite.fr},
#'         Nicolas Loiseau, \email{nicolas.loiseau1@@gmail.com},
#'         Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr}
#'
#' @date 2020/06/11


owarn <- options()$warn
options(warn = -1)



for (horizon in cc_horizons) {



#'  -------------------------------------------------------------------------   @InitExportDevice


  png(
    file      = here::here("figures", paste0(figname, "_", horizon,".png")),
    width     = 24.00,
    height    = 12.00,
    units     = "in",
    res       = 600,
    pointsize = 38
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
    cex.axis = 0.70
  )

  par(mfcol = c(length(taxas), length(threats_vars) + 1))


  comp <- 1

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


    grid(lty = 1, lwd = .5)
    box(lwd = 2)



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
      axis(side = 1, lwd = 0, font = 2, cex.axis = 0.75)

    }



#'  -------------------------------------------------------------------------   @AddyAxis


    par(mgp = c(0.5, 0, 0))
    axis(side = 2, at = 1:3, labels = paste(classes_labs, c("†", "¶", "§")), lwd = 0, font = 2, cex.axis = 0.85)



#'  -------------------------------------------------------------------------   @AddTitle


    if (taxa == "mammals") {

      mtext(text = "IUCN status (%)", side = 3, font = 2, cex = 0.75)

    }



#'  -------------------------------------------------------------------------   @AddSilhouette


    addPhylopic(
      img   = icons[[taxa]],
      x     = 0.11,
      y     = 0.93,
      ysize = 0.20,
      alpha = 1,
      color = color_silh,
      AR    = 0.6
    )

    text(
      x      = par()$usr[2] - 0.05 * (par()$usr[2] - par()$usr[1]),
      y      = par()$usr[4] - 0.04 * (par()$usr[4] - par()$usr[3]),
      labels = letters[comp],
      col    = dark_grey,
      font   = 2,
      cex    = 1
    )

    comp <- comp + 1

  } # e_o taxas



#'  -------------------------------------------------------------------------   @ThreatsPlots



  for (threat in threats_vars) {



    if (threat == "climate_change") {

      colname <- paste0("cc_", horizon)

    } else {

      colname <- threat

    }

    subdatas <- datas[!is.na(datas[ , colname]), c("class", "dr_class", colname)]
    subdatas <- subdatas[subdatas[ , "dr_class"] %in% classes, ]

    if (threat == "targetmet_percentagecover") {

      subdatas <- subdatas[subdatas[ , threat] <= 100, ]

    }

    for (taxa in taxas) {



      if (comp / 2 == round(comp / 2)) {

        par(mar = c(1.0, 1.0, 0.5, 0.5))

      } else {

        par(mar = c(0.5, 1.0, 1.0, 0.5))

      }

      ssubdatas <- subdatas[subdatas[ , "class"] == taxa, ]



#'  -------------------------------------------------------------------------   @ComputePvalue

      dat     <- ssubdatas
      colnames(dat)[ncol(dat)] <- "values"

      mod     <- aov(values ~ dr_class, data = dat)
      posthoc <- TukeyHSD(x = mod, which = "dr_class", conf.level = 0.95)
      posthoc <- posthoc$dr_class

      vars <- strsplit(rownames(posthoc), "-")

      pval <- data.frame(
        taxa             = taxa,
        threat           = colname,
        var1             = unlist(lapply(vars, function(x) x[1])),
        var2             = unlist(lapply(vars, function(x) x[2])),
        pval             = posthoc[ , "p adj"],
        row.names        = NULL,
        stringsAsFactors = FALSE
      )

      dat <- data.frame()
      for (i in 1:nrow(pval)) {

        tmp <- rbind(pval[i, ], pval[i, ])
        dat <- rbind(dat, tmp)

      }

      for (i in seq(2, nrow(dat), 2)) {

        xx <- dat[i, "var1"] ; yy <- dat[i, "var2"]
        dat[i, "var1"] <- yy ; dat[i, "var2"] <- xx

      }

      pval <- data.frame()

      for (classe in classes) {

        tmp <- dat[dat[ , "var1"] == classe , ]

        xxx <- NULL

        for (i in 1:nrow(tmp)) {

          if (tmp[i, "pval"] >= 0.05) { xxx <- c(xxx, tmp[i, "var2"]) }

        }

        xxx <- gsub("D75R75", "§", xxx)
        xxx <- gsub("AVG", "¶", xxx)
        xxx <- gsub("D25R25", "†", xxx)
        # xxx <- gsub("D75R75", "*", xxx)
        # xxx <- gsub("AVG", "*", xxx)
        # xxx <- gsub("D25R25", "*", xxx)

        xxx <- paste0(sort(xxx), collapse = " ")

        tmp <- data.frame(
          taxa             = taxa,
          threat           = colname,
          dr_class         = classe,
          signif           = xxx,
          row.names        = NULL,
          stringsAsFactors = FALSE
        )

        pval <- rbind(pval, tmp)
      }


#'  -------------------------------------------------------------------------   @EmptyPlot


      xmin <- min(subdatas[ , colname], na.rm = TRUE)
      xmax <- max(subdatas[ , colname], na.rm = TRUE)

      xmax <- xmax + (xmax - xmin) * 0.22

      plot(0,
        xlim = c(xmin, xmax),
        ylim = c(0.5, 3.5),
        axes = FALSE,
        ann  = FALSE,
        type = "n"
      )



#'  -------------------------------------------------------------------------   @AddBackground


      grid(lty = 1, lwd = .5)
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


        # boxplot(
        #   x          = x_coords,
        # 	names      = NULL,
        # 	horizontal = TRUE,
        # 	col        = "transparent",
        #   border     = "black",
        # 	lty        = 1,
        # 	lwd        = 1,
        # 	at         = which(classes == classe),
        # 	add        = TRUE,
        #   outline    = FALSE,
        #   axes       = FALSE,
        #   width      = 1.8
        # )
        vioplot::vioplot(
          x           = x_coords,
          names       = NULL,
          horizontal  = TRUE,
          col         = "transparent",
          border      = dark_grey,
          lty         = 1,
          lwd         = 2,
          rectCol     = dark_grey,
          lineCol     = dark_grey,
          colMed      = "white",
          colMed2     = "white",
          at          = which(classes == classe),
          add         = TRUE,
          na.rm       = TRUE,
          side        = "both",
          plotCentre  = "points"
        )


#'  -------------------------------------------------------------------------   @AddVerticalLineCC


        if (threat == "climate_change") {

          abline(v = 0, lwd = 1, lty = 2)

        }



#'  -------------------------------------------------------------------------   @AddPvalues


        p_val <- pval[which(pval[ , "dr_class"] == classe), "signif"]

        if (length(pval) > 0) {

          x_coords <- ssubdatas[ , colname]
          x_coords <- x_coords[!is.na(x_coords)]

          text(
            x      = max(x_coords),
            y      = which(classes == classe),
            labels = p_val,
            cex    = 1.00,
            font   = 2,
            pos    = 4,
            col = dark_grey
          )

        }

      } # e_o classes



#'  -------------------------------------------------------------------------   @AddxAxis


      if (comp / 2 == round(comp / 2)) {

        par(mgp = c(0.5, -0.05, 0))
        axis(side = 1, lwd = 0, font = 2, cex.axis = 0.75)

      }



#'  -------------------------------------------------------------------------   @AddTitle


      if (comp / 2 != round(comp / 2)) {

        mtext(
          text = threats_labs[which(threats_vars == threat)],
          side = 3,
          font = 2,
          cex  = 0.75
        )

      }

      text(
        x      = par()$usr[2] - 0.05 * (par()$usr[2] - par()$usr[1]),
        y      = par()$usr[4] - 0.04 * (par()$usr[4] - par()$usr[3]),
        labels = letters[comp],
        col    = dark_grey,
        font   = 2,
        cex    = 1
      )

      comp <- comp + 1

    } # e_o taxas
  } # e_o threats_vars


  dev.off()

  ## Message ----

  usethis::ui_done(
    paste(
      usethis::ui_value(figname),
      paste0("(", gsub("_", "-", horizon), ")"),
      "successfully exported in figures/"
    )
  )
} # e_o cc_horizons

options(warn = owarn)
