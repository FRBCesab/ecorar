#' --------------------------------------------------------------------------   @Header
#'
#' @title Run project
#'
#' @description
#' Run project.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @date 2019/10/25
#'
#' --------------------------------------------------------------------------   @Header



rm(list = ls())



#'  -------------------------------------------------------------------------   @ProjectSetup


source(file.path("analyses", "setup.R"))



#'  -------------------------------------------------------------------------   @Figure1


source(file.path("analyses", "figure1_pcoa.R"))



#'  -------------------------------------------------------------------------   @Figure2


source(file.path("analyses", "figure2_phylo.R"))



#'  -------------------------------------------------------------------------   @Figure3


source(file.path("analyses", "figure3_maps.R"))



#'  -------------------------------------------------------------------------   @Figure4


source(file.path("analyses", "figure4_boxplots.R"))
