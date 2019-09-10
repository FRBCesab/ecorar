################################################################################
###                                                                          ###
###                      FIGURE 4 (THREATS VIOLIN PLOTS)                     ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### AUTHORS : The Three Nico                                                 ###
### DATE    : 2019/08/05                                                     ###
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

filename <- "Figure_4"

horizon  <- "2061-2080"



#'  -------------------------------------------------------------------------   @LoadAddings


library(vioplot)
library(png)

addings <- list.files(path = path_R, pattern = "\\.R$", full.names = TRUE)
for (i in 1:length(addings)) { source(addings[i]) }




#'  -------------------------------------------------------------------------   @ImportSpeciesSilhouettes


icons <- list()

for (taxa in taxas) {

  icons[[taxa]] <- readPNG(
      source = file.path(path_data, taxa, paste0(taxa, "_silhouette.png"))
  )
}



#'  -------------------------------------------------------------------------   @ImportData


threats <- get(
  load(
    file = file.path(path_data, "threats_nocc.RData")
  )
)

for (i in 2:ncol(threats)) { threats[ , i] <- as.character(threats[ , i]) }


climate <- get(
  load(
    file = file.path(path_data, "threats_cc.RData")
  )
)

for (i in 2:ncol(climate)) { climate[ , i] <- as.character(climate[ , i]) }

climate <- climate[climate[ , "Horizon"] == horizon, -ncol(climate)]


iucn <- get(
  load(
    file = file.path(path_data, "iucn_status.RData")
  )
)

pvalues <- get(
  load(
    file = file.path(path_data, paste0("violin_pvalues_", horizon, ".RData"))
  )
)


#'  -------------------------------------------------------------------------   @PrepareData


datas <- rbind(threats, climate)

values_m  <- tapply(
  X      = datas$Value,
  INDEX  = list(datas$Threats, datas$Taxa, datas$DR_class),
  FUN    = function(x) { mean(x, na.rm = TRUE) }
)

values_sd <- tapply(
  X      = datas$Value,
  INDEX  = list(datas$Threats, datas$Taxa, datas$DR_class),
  FUN    = function(x) { sd(x, na.rm = TRUE) }
)

iucn <- lapply(
  X    = iucn,
  FUN  = function(y) { apply(y, 2, function(x) x / sum(x) * 100) }
)



#'  -------------------------------------------------------------------------   @InitExportDevice


png(
  file       = file.path(path_figs, paste0(filename, "_", horizon, ".png")),
  width      = 12.00,
  height     =  7.00,
  units      = "in",
  res        = 600,
  pointsize  = 18
)



#'  -------------------------------------------------------------------------   @DefineGraphicalParameters


par(
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
  cex.axis  = 0.75,
  fg        = "#666666"
)

par(mfcol = c(2, 5))



for (taxa in taxas) {



  if (taxa == "birds") {   par(mar = c(1.0, 1.0, 0.5, 0.5)) }

  if (taxa == "mammals") { par(mar = c(0.5, 1.0, 1.0, 0.5)) }



#'  -------------------------------------------------------------------------   @EmptyPlotIUCN


  plot(0,
    xlim  = c(-10, 110),
    ylim  = c(0.5, 3.5),
    axes  = FALSE,
    type  = "n",
    ann   = FALSE
  )



#'  -------------------------------------------------------------------------   @AddBackground


  grid()
  box()



#'  -------------------------------------------------------------------------   @AddIUCNStickers


  for (classe in classes) {

    for (status in iucn_status) {

      points(
        x    = iucn[[taxa]][status, classe],
        y    = which(classes == classe),
        cex  = 4,
        col  = color_iucn_bg[status],
        pch  = 19
      )

      text(
        x       = iucn[[taxa]][status, classe],
        y       = which(classes == classe),
        labels  = status,
        cex     = 0.85,
        col     = color_iucn_fg[status],
        family  = "sans",
        font    = 2
      )
    }
  }



#'  -------------------------------------------------------------------------   @AddxAxis


  if (taxa == "birds") {

    par(mgp = c(0.5, -0.05, 0))
    axis(1, lwd = 0, font = 2)
  }



#'  -------------------------------------------------------------------------   @AddyAxis


    par(mgp = c(0.5, 0, 0))
    axis(2, at = 1:3, labels = classes_labs, lwd = 0, font = 2)



#'  -------------------------------------------------------------------------   @AddTitle


  if (taxa == "mammals") {

    mtext(text = "IUCN status (%)", side = 3, font = 2, cex = 0.85)

  }



#'  -------------------------------------------------------------------------   @AddSilhouette


  addPhylopic(
    img    = icons[[taxa]],
    x      = 0.85,
    y      = 0.92,
    ysize  = 0.20,
    alpha  = 1,
    color  = "#777777",
    AR     = 0.5
  )
}



#'  -------------------------------------------------------------------------   @ThreatsPlots


comp <- 3



for (threat in threats_vars) {



  tmp0 <- datas[which(datas[ , "Threats"] == threat), ]



  for (taxa in taxas) {



    if (comp / 2 == round(comp / 2)) {

      par(mar = c(1.0, 1.0, 0.5, 0.5))

    } else {

      par(mar = c(0.5, 1.0, 1.0, 0.5))
    }

    tmp <- tmp0[which(tmp0[ , "Taxa"] == taxa), ]



#'  -------------------------------------------------------------------------   @EmptyPlot


    xmin <- min(tmp0[ , "Value"])
    xmax <- max(tmp0[ , "Value"])
    xmax <- xmax + (xmax - xmin) * 0.15

    plot(0,
      xlim  = c(xmin, xmax),
      ylim  = c(0.5, 3.5),
      axes  = FALSE,
      type  = "n",
      ann   = FALSE
    )



#'  -------------------------------------------------------------------------   @AddBackground


    grid()
    box()

    if (comp > 8) { abline(v = 0, lwd = 1, lty = 2) }



    for (classe in classes){



#'  -------------------------------------------------------------------------   @AddViolonPlot


      vioplot(
        x           = tmp[which(tmp[ , "DR_class"] == classe), "Value"],
      	names       = NULL,
      	horizontal  = TRUE,
      	col         = paste0(color_classes[classe], "BB"),
        border      = "#666666",
      	lty         = 1,
      	lwd         = 1,
      	rectCol     = color_classes[classe],
        lineCol     = color_classes[classe],
      	colMed      = color_classes[classe],
        colMed2     = "white",
      	at          = which(classes == classe),
      	add         = TRUE,
      	na.rm       = TRUE,
      	side        = "both",
        plotCentre  = "line"
      )



#'  -------------------------------------------------------------------------   @AddViolonPlot


      val <- pvalues[
        which(
          pvalues[ , "taxa"] == taxa &
          pvalues[ , "threat"] == threat &
          pvalues[ , "DR_class"] == classe
        ), "signif"]

      text(
        x       = max(tmp[which(tmp[ , "DR_class"] == classe), "Value"]),
        y       = which(classes == classe),
        labels  = gsub(",", ", ", val),
        cex     = 0.85,
        font    = 2,
        pos     = 4
      )
    }



#'  -------------------------------------------------------------------------   @AddxAxis


    if (comp / 2 == round(comp / 2)) {

      par(mgp = c(0.5, -0.05, 0))
      axis(1, lwd = 0, font = 2)
    }



#'  -------------------------------------------------------------------------   @AddTitle


    if (comp / 2 != round(comp / 2)) {

      mtext(
        text  = threats_labels[which(threats_vars == threat)],
        side  = 3,
        font  = 2,
        cex   = 0.85
      )
    }

    comp <- comp + 1
  }
}

dev.off()
