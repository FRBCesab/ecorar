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



### PATH TO DATA
path_data <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/data/"       ### !!!


### IMPORT DATA
birds   <- readRDS(paste0(path_data, "funk_birds.rds"))
mammals <- readRDS(paste0(path_data, "funk_mammals.rds"))
