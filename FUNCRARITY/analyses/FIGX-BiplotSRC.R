################################################################################
###                                                                          ###
###                       FIGURE X (BIPLOT SRC vs AREA)                      ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### AUTHORS : The Three Nico                                                 ###
### DATE    : 2019/08/07                                                     ###
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


root      <- "/Users/nicolascasajus/OneDrive/OneDrive - Fondation Biodiversité/MySpace/GROUPS/FREE/01-Loiseau/FUNCRARITY/"
source(file.path(root, "graphsParameters.R"))

filename <- "Figure_X"

horizon  <- "2080"



#'  -------------------------------------------------------------------------   @LoadAddings


addings <- list.files(path = path_R, pattern = "\\.R$", full.names = TRUE)
for (i in 1:length(addings)) { source(addings[i]) }



#'  -------------------------------------------------------------------------   @ImportData


src <- get(
  load(
    file = file.path(path_data, "src_vs_occurrence.RData")
  )
)

src[ , "occurrence_log"] <- log(src[ , "occurrence"])
src[ , "horizon"]        <- src[ , paste0("src", horizon)]



#'  -------------------------------------------------------------------------   @InitExportDevice


png(
  file       = file.path(path_figs, paste0(filename, "_", horizon, ".png")),
  width      = 12.00,
  height     = 12.00,
  units      = "in",
  res        = 600,
  pointsize  = 18
)



#'  -------------------------------------------------------------------------   @DefineGraphicalParameters


par(
  family    = "serif",
  mar       = c(2.75, 2.75, 0.5, 0.5),
  cex.axis  = 1.25,
  mgp       = c(2, .25, 0),
  tcl       = -0.25,
  xpd       = FALSE,
  new       = FALSE,
  fig       = c(0, 1, 0, 1),
  col       = "#666666",
  col.axis  = "#666666",
  cex.axis  = 0.75,
  fg        = "#666666"
)



#'  -------------------------------------------------------------------------   @EmptyPlot


plot(0,
  xlim  = range(src[ , "occurrence_log"]),
  ylim  = range(src[ , paste0("src", horizon)]),
  axes  = FALSE,
  type  = "n",
  ann   = FALSE
)

grid()

points(
  x    = src[src[ , "taxa"] == "birds", "occurrence_log"],
  y    = src[src[ , "taxa"] == "birds", paste0("src", horizon)],
  pch  = 19,
  col  = paste0(color_avg, "88"),
  cex  = 1
)

points(
  x    = src[src[ , "taxa"] == "mammals", "occurrence_log"],
  y    = src[src[ , "taxa"] == "mammals", paste0("src", horizon)],
  pch  = 19,
  col  = paste0(color_rare, "88"),
  cex  = 1
)


xpred <- data.frame(
  occurrence_log = seq(
    min(src[ , "occurrence_log"]),
    max(src[ , "occurrence_log"]),
    length.out = 1000
  )
)

src_birds   <- src[src[ , "taxa"] == "birds", ]
src_mammals <- src[src[ , "taxa"] == "mammals", ]

mod_birds   <- lm(horizon ~ occurrence_log, data = src_birds)
mod_mammals <- lm(horizon ~ occurrence_log, data = src_mammals)

y_birds   <- predict(mod_birds, xpred, type = "response", se.fit = TRUE)
y_mammals <- predict(mod_mammals, xpred, type = "response", se.fit = TRUE)


polygon(
  x      = c(xpred$occurrence_log, xpred$occurrence_log[length(xpred$occurrence_log):1]),
  y      = c(y_birds$fit + y_birds$se.fit, (y_birds$fit - y_birds$se.fit)[length(y_birds$fit):1]),
  col    = "#cdcdcdaa",
  border = "#cdcdcdaa"
)

# Add smooth prediction
lines(
  x   = xpred$occurrence_log,
  y   = y_birds$fit,
  col = color_avg,
  lwd = 2.25,
  lty = 1
)

polygon(
  x      = c(xpred$occurrence_log, xpred$occurrence_log[length(xpred$occurrence_log):1]),
  y      = c(y_mammals$fit + y_mammals$se.fit, (y_mammals$fit - y_mammals$se.fit)[length(y_mammals$fit):1]),
  col    = "#cdcdcdaa",
  border = "#cdcdcdaa"
)

# Add smooth prediction
lines(
  x   = xpred$occurrence_log,
  y   = y_mammals$fit,
  col = color_rare,
  lwd = 2.25,
  lty = 1
)



axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
mtext(side = 1, line = 1.5, "Number of cells (log)", font = 2)
mtext(side = 2, line = 1.5, paste0("Species range change (CC ", horizon, ")"), font = 2)
box()

abline(h = 0, lty = 2)


rect(8.1, 455, 8.6, 465,
  col    = "#cdcdcdaa",
  border = "#cdcdcdaa"
)

lines(
  x = c(8.1, 8.6),
  y = rep(460, 2),
  col = color_avg,
  lwd = 2.25,
  lty = 1
)

rect(8.1, 435, 8.6, 445,
  col    = "#cdcdcdaa",
  border = "#cdcdcdaa"
)

lines(
  x = c(8.1, 8.6),
  y = rep(440, 2),
  col = color_rare,
  lwd = 2.25,
  lty = 1
)

text(
  x = 8.6,
  y = 460,
  labels = "Birds",
  font = 2,
  pos = 4
)

text(
  x = 8.6,
  y = 440,
  labels = "Mammals",
  font = 2,
  pos = 4
)


dev.off()
