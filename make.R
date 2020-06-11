#' Run the Entire Project
#'
#' This script runs the entire project and produces all figures presents in the
#' Loiseau, Mouquet et al.'s 2020 paper.
#'
#' @author Nicolas Casajus, \email{nicolas.casajus@@fondationbiodiversite.fr},
#'         Nicolas Loiseau, \email{nicolas.loiseau1@@gmail.com},
#'         Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr}
#'
#' @date 2020/06/10


if (!("here" %in% installed.packages())) install.packages("here")

source(here::here("analyses", "setup.R"))
source(here::here("analyses", "params.R"))
source(here::here("analyses", "imports.R"))

# figname <- "Figure_1"
# source(here::here("analyses", "figure_biplot.R"))
# 
# figname <- "Figure_2"
# source(here::here("analyses", "figure_pcoa.R"))
# 
# figname <- "Figure_3"
# source(here::here("analyses", "figure_phylo.R"))

figname <- "Figure_4"
source(here::here("analyses", "figure_maps.R"))

