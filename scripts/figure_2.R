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
### [1] rgdal_1.4-3   rgeos_0.4-2   raster_2.8-19    sp_1.3-1                ###
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



### PATH TO DATA
path_data <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/"       ### !!!
path_R    <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/RALLL/"       ### !!!

source(paste0(path_R, "scripts/__addGraticules.R"))
source(paste0(path_R, "scripts/__compassRose.R"))

### IMPORT DATA
study   <- raster(paste0(path_data, "data/reference_grid_50km.tif"))
birds   <- readRDS(paste0(path_data, "data/funk_birds.rds"))
mammals <- readRDS(paste0(path_data, "data/funk_mammals.rds"))
world   <- readOGR(dsn = paste0(path_data, "data/ne_10m_ocean"), layer = "ne_10m_ocean")
world   <- spTransform(world, proj4string(study))

### DEFINE MAP EXTENT (WITH ANTARCTICA)
border <- as(extent(c(-180, 180, -90, 90)), "SpatialPolygons")
proj4string(border) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
border <- spTransform(border, proj4string(birds))


### CREATE OCEAN LAYER
ocean <- study
ocean[][which(!is.na(ocean[]))] <- 0
ocean[][which(is.na(ocean[]))] <- 1
ocean[][which(ocean[] == 0)] <- NA


### COLOR GRADIENT
cols <- c("#aaaaaa", colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(255))


### SPACING BOX (INSET)
padding <- 25000


### SINGLE MAP
png(
  file   = paste0("~/Desktop/test-mammals.png"),
  width  = 12,
  height = 5.25,
  units  = "in",
  res    = 600
)

par(
  xaxs      = "i",
  yaxs      = "i",
  family    = "serif",
  mar       = rep(0.75, 4),
  cex.axis  = 1.25,
  mgp       = c(2, .5, 0),
  tcl       = -.25,
  xpd       = FALSE,
  new       = FALSE,
  fig       = c(0, 1, 0, 1),
  col       = "#666666",
  col.axis  = "#666666",
  fg        = "#666666"
)

plot(border, border = NA, col = "#95D8EB")

plot(subset(mammals, "mammals_D75R75"), col = cols, legend = FALSE, add = TRUE)

plot(world, col = "#95D8EB", border = NA, add = TRUE)


zzz = addGraticules(add = TRUE, line.color = "#aaaaaa")

par(mgp = c(.5, -.1, 0))
axis(1, at = zzz[[1]][ , "x"], labels = paste0(zzz[[1]][ , "label"], "°", zzz[[1]][ , "direction"]), cex.axis = .75, lwd = 0)
par(mgp = c(.5, 0, 0))
axis(3, at = zzz[[1]][ , "x"], labels = paste0(zzz[[1]][ , "label"], "°", zzz[[1]][ , "direction"]), cex.axis = .75, lwd = 0)

par(mgp = c(.5, 0, 0))
axis(2, at = zzz[[2]][ , "y"], labels = paste0(zzz[[2]][ , "label"], "°", zzz[[2]][ , "direction"]), cex.axis = .75, lwd = 0)
par(mgp = c(.5, -.1, 0))
axis(4, at = zzz[[2]][ , "y"], labels = paste0(zzz[[2]][ , "label"], "°", zzz[[2]][ , "direction"]), cex.axis = .75, lwd = 0)


par(xpd = TRUE)
plot(border, border = par()$col, lwd = 4, add = TRUE)
plot(border, border = "white", lwd = 2, add = TRUE)
par(xpd = FALSE)

cadre1 <- c(xmn = -8879810, xmx = -5516986, ymn = 632368.8, ymx = 3118863.9)
rect(cadre1["xmn"], cadre1["ymn"], cadre1["xmx"], cadre1["ymx"], border = "#888888", lwd = 1)

lines(
  x    = rep(par()$usr[1] + par()$usr[2] * 0.25, 2),
  y    = c(cadre1["ymn"] + (cadre1["ymx"] - cadre1["ymn"]) / 2, par()$usr[3] + par()$usr[4] * 0.9475),
  col  = "#888888"
)

lines(
  x    = c(cadre1["xmn"], par()$usr[1] + par()$usr[2] * 0.25),
  y    = rep(cadre1["ymn"] + (cadre1["ymx"] - cadre1["ymn"]) / 2, 2),
  col  = "#888888"
)


cadre2 <- c(xmn = 3980550, xmx = 5215893, ymn = -997417, ymx = -3517876)
rect(cadre2["xmn"], cadre2["ymn"], cadre2["xmx"], cadre2["ymx"], border = "#888888", lwd = 1)

lines(
  x    = rep(cadre2["xmn"] + (cadre2["xmx"] - cadre2["xmn"]) / 2, 2),
  y    = c(cadre2["ymx"], par()$usr[3] + par()$usr[4] * 0.25),
  col  = "#888888"
)

lines(
  x    = c(cadre2["xmn"] + (cadre2["xmx"] - cadre2["xmn"]) / 2, par()$usr[2] * 0.01),
  y    = rep(par()$usr[3] + par()$usr[4] * 0.25, 2),
  col  = "#888888"
)


cadre3 <- c(xmn = 9047080, xmx = 15780354, ymn = -1573001, ymx = 903988)
rect(cadre3["xmn"], cadre3["ymn"], cadre3["xmx"], cadre3["ymx"], border = "#888888", lwd = 1)

lines(
  x    = rep(par()$usr[2] * 0.475, 2),
  y    = c(cadre3["ymn"] + (cadre3["ymx"] - cadre3["ymn"]) / 2, -1 * par()$usr[4] * 0.5),
  col  = "#888888"
)

lines(
  x    = c(cadre3["xmn"], par()$usr[2] * 0.475),
  y    = rep(cadre3["ymn"] + (cadre3["ymx"] - cadre3["ymn"]) / 2, 2),
  col  = "#888888"
)


cadre4 <- c(xmn = 11053993, xmx = 12457308, ymn = 420761, ymx = 3458091)
rect(cadre4["xmn"], cadre4["ymn"], cadre4["xmx"], cadre4["ymx"], border = "#888888", lwd = 1)
#
# lines(
#   x    = rep(par()$usr[2] * 0.475, 2),
#   y    = c(cadre4["ymn"] + (cadre4["ymx"] - cadre4["ymn"]) / 2, -1 * par()$usr[4] * 0.5),
#   col  = "#888888"
# )
#
lines(
  x    = c(cadre4["xmx"], par()$usr[2] * 0.90),
  y    = rep(cadre4["ymn"] + (cadre4["ymx"] - cadre4["ymn"]) / 2, 2),
  col  = "#888888"
)

xstart <- xleft <- -17050000
ybottom <-  -6550000

for (i in 1:length(cols)) {
  xleft   <- xleft + 35000
  rect(xleft, ybottom, xleft + 35000, ybottom + 350000, border = NA, col = cols[i])
}
rect(xstart, ybottom, xleft + 35000, ybottom + 350000, border = "white", col = NA)
text(xstart, ybottom, "0", pos = 1, col = par()$col.axis, font = 2)
text(xleft + 35000, ybottom, max(subset(mammals, "mammals_D75R75")[], na.rm = TRUE), pos = 1, col = par()$col.axis, font = 2)
text(xstart + (((xleft + 35000) - xstart) / 2), ybottom + 350000, "Number of mammal species (D75R75)", pos = 3, col = par()$col.axis, font = 2)
text(xstart + (((xleft + 35000) - xstart) / 2), ybottom, max(subset(mammals, "mammals_D75R75")[], na.rm = TRUE) / 2, pos = 1, col = par()$col.axis, font = 2)


compassRose(-14000000, 4900000, cex.cr = .65, cex.let = .65)




par(new = TRUE, fig = c(0.05, 0.25, 0.15, 0.50))

plot(0, xlim = c(cadre1["xmn"], cadre1["xmx"]), ylim = c(cadre1["ymn"], cadre1["ymx"]), type = "n", axes = FALSE, bty = "n", ann = FALSE)

par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = "#c0e8f3",
  border   = par()$col,
  lwd      = 1
)
par(xpd = FALSE)


plot(subset(mammals, "mammals_D75R75"), col = cols, legend = FALSE, axes = FALSE, box = FALSE, add = TRUE)


par(fig = c(0.05, 0.25, 0.15, 0.50))

par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = par()$col,
  lwd      = 3
)
par(xpd = FALSE)


plot(world, col = "#c0e8f3", border = NA, add = TRUE)
zzz = addGraticules(add = TRUE)
par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = par()$col,
  lwd      = 3
)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = "white",
  lwd      = 1.5
)
par(xpd = FALSE)
text(par()$usr[2] - 225000, par()$usr[4] - 200000, "a)", col = par()$col.axis, font = 2)

par(new = TRUE)
par(fig = c(0.41, 0.52, 0.025, 0.490))

plot(0, xlim = range(c(cadre2["xmn"], cadre2["xmx"])), ylim = range(c(cadre2["ymn"], cadre2["ymx"])), type = "n", axes = FALSE, bty = "n", ann = FALSE)

par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = "#c0e8f3",
  border   = par()$col,
  lwd      = 1
)
par(xpd = FALSE)

plot(subset(mammals, "mammals_D75R75"), col = cols, legend = FALSE, axes = FALSE, box = FALSE, add = TRUE)


par(fig = c(0.41, 0.52, 0.025, 0.490))

par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = par()$col,
  lwd      = 3
)
par(xpd = FALSE)


plot(world, col = "#c0e8f3", border = NA, add = TRUE)
zzz = addGraticules(add = TRUE)
par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = par()$col,
  lwd      = 3
)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = "white",
  lwd      = 1.5
)
par(xpd = FALSE)
text(par()$usr[2] - 150000, par()$usr[4] - 150000, "b)", col = par()$col.axis, font = 2)

par(new = TRUE)
par(fig = c(0.67, 0.98, 0.025, 0.325))

plot(0, xlim = range(c(cadre3["xmn"], cadre3["xmx"])), ylim = range(c(cadre3["ymn"], cadre3["ymx"])), type = "n", axes = FALSE, bty = "n", ann = FALSE)

par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = "#c0e8f3",
  border   = par()$col,
  lwd      = 1
)
par(xpd = FALSE)

plot(subset(mammals, "mammals_D75R75"), col = cols, legend = FALSE, axes = FALSE, box = FALSE, add = TRUE)

par(fig = c(0.67, 0.98, 0.025, 0.325))

par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = par()$col,
  lwd      = 3
)
par(xpd = FALSE)


plot(world, col = "#c0e8f3", border = NA, add = TRUE)
zzz = addGraticules(add = TRUE)
par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = par()$col,
  lwd      = 3
)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = "white",
  lwd      = 1.5
)
par(xpd = FALSE)
text(par()$usr[2] - 250000, par()$usr[4] - 300000, "c)", col = par()$col.axis, font = 2)


par(new = TRUE)
par(fig = c(0.88, 0.98, 0.475, 0.97))

plot(0, xlim = range(c(cadre4["xmn"], cadre4["xmx"])), ylim = range(c(cadre4["ymn"], cadre4["ymx"])), type = "n", axes = FALSE, bty = "n", ann = FALSE)

par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = "#c0e8f3",
  border   = par()$col,
  lwd      = 1
)
par(xpd = FALSE)

plot(subset(mammals, "mammals_D75R75"), col = cols, legend = FALSE, axes = FALSE, box = FALSE, add = TRUE)

par(fig = c(0.88, 0.98, 0.475, 0.97))

par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = par()$col,
  lwd      = 3
)
par(xpd = FALSE)


plot(world, col = "#c0e8f3", border = NA, add = TRUE)
zzz = addGraticules(add = TRUE)
par(xpd = TRUE)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = par()$col,
  lwd      = 3
)
rect(
  xleft    = par()$usr[1] - padding,
  ybottom  = par()$usr[3] - padding,
  xright   = par()$usr[2] + padding,
  ytop     = par()$usr[4] + padding,
  col      = NA,
  border   = "white",
  lwd      = 1.5
)

text(par()$usr[2] - 200000, par()$usr[4] - 200000, "d)", col = par()$col.axis, font = 2)
par(xpd = FALSE)




dev.off()
