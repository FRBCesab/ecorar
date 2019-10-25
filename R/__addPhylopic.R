#' @title Add transparency and/or change color of a picture
#'
#' @description
#' This function recolours and plots a picture. Not perfect, needs improvement.
#'
#' @param img [array] a picture loaded by \code{png::readPNG()}.
#' @param x [numeric] the x coordinate of the center of the picture.
#' @param y [numeric] the y coordinate of the center of the picture.
#' @param ysize [numeric] a facteur to resize picture dimensions.
#' @param alpha [numeric] transparency factor between 0 (transparent) and 1 (opaque).
#' @param color [array] a color name (e.g. "red", "steelblue") or an hexadecimal value (e.g. "#ff0000")
#' @param AR [array] a magnification factor (zoom).
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#'
#' @export
#'
#' @return
#' This function returns a picture matrix with colors and/or transparency altered.
#' This matrix can be plotted with the function \code{graphics::rasterImage()} (see examples).
#'
#' @examples
#'
#' # No example.



addPhylopic <- function(
  img, x = NULL, y = NULL, ysize = NULL, alpha = 0.2, color = NULL, AR = 1
) {


  img <- recolorPhylopic(img, alpha, color)

  graphics::par(usr = c(0, 1, 0, 1))

  graphics::rasterImage(
    image        = img,
    xleft        = x - (ysize / 2),
    ybottom      = y - (AR * ysize / 2),
    xright       = x + (ysize / 2),
    ytop         = y + (AR * ysize / 2),
    interpolate  = TRUE
  )
}
