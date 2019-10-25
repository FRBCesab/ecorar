#' --------------------------------------------------------------------------   @Header
#'
#' @title Project setup
#'
#' @description
#' This R script...
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @date 2019/10/25
#'
#' --------------------------------------------------------------------------   @Header



rm(list = ls())



#' ----------------------------------------------------------------------------- @InstallCranLibs


pkgs <- c(
  "png",
  "ggplot2",
  "cowplot",
  "grid",
  "gridExtra",
  "sp",
  "rgdal",
  "rgeos",
  "raster",
  "RColorBrewer"
)

nip <- pkgs[!(pkgs %in% installed.packages())]

nip <- lapply(nip, install.packages, dependencies = TRUE)



#' ----------------------------------------------------------------------------- @InstallDevLibs



#' ----------------------------------------------------------------------------- @LoadLibs


ip <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

if (sum(ip) != length(pkgs)) { cat("Some packages failed to load.\n") }



#' ----------------------------------------------------------------------------- @LoadRFunctions

rfun <- list.files(path = "R", pattern = "^__.+\\.R$", full.names = TRUE)

rfun <- unlist(lapply(rfun, source, verbose = FALSE))



#' ----------------------------------------------------------------------------- @CreateFolders


dir_names <- "figures"
dir_vars  <- "path_figs"

dirs <- lapply(

  X   = 1:length(dir_names),

  FUN = function(i) {

    dir.create(
      path         = dir_names[i],
      showWarnings = FALSE,
      recursive    = TRUE
    )

    assign(
      x     = dir_vars[i],
      value = dir_names[i],
      envir = .GlobalEnv
    )
  }
)


rm(list = c("dirs", "dir_names", "dir_vars", "pkgs", "nip", "ip", "rfun"))



#' ---------------------------------------------------------------------------- @GeneralParameters
