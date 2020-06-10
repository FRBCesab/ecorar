#' @title Get corners of convex hull of a set of points
#'
#' @description
#' This function gets corners coordinates of the convex hull of a set of points.
#'
#' @param data [data.frame] a data frame with at least two columns named x and y.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#' @author Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @export
#'
#' @return
#' This function returns a susbet of the initial data frame where each row is
#' a corner of the convex hull.
#'
#' @examples
#'
#' x   <- rnorm(100)
#' y   <- rnorm(100)
#' mat <- data.frame(x, y)
#'
#' dat <- findHull(mat)
#'
#' dat <- rbind(dat, dat[1, ])
#'
#' plot(mat, pch = 19)
#' polygon(dat[ , 1], dat[ , 2], border = "red")



findHull <- function(data) {


  if (sum(which(colnames(data) == "x")) == 0) { stop("No column `x` found in `data`.") }
  if (sum(which(colnames(data) == "y")) == 0) { stop("No column `y` found in `data`.") }

  if (!is.numeric(data$x)) { stop("`x` must be numeric`.") }
  if (!is.numeric(data$y)) { stop("`y` must be numeric`.") }

  return(data[chull(data$x, data$y), ])

}
