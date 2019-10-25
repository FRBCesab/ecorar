#' --------------------------------------------------------------------------   @Header
#'
#' @title Run project
#'
#' @description
#' Run project.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#'
#' @date 2019/10/25
#'
#' --------------------------------------------------------------------------   @Header



rm(list = ls())



#'  -------------------------------------------------------------------------   @ProjectSetup


source(file.path("analysis", "setup.R"))



#'  -------------------------------------------------------------------------   @Figure1


source(file.path("analysis", "figure1_pcoa.R"))



#'  -------------------------------------------------------------------------   @Figure2


source(file.path("analysis", "figure2_pcoa.R"))



#'  -------------------------------------------------------------------------   @Figure3


source(file.path("analysis", "figure3_pcoa.R"))
