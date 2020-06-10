#' @title Add transparency and/or change color of a picture
#'
#' @description
#' This function adds transparency and/or change color of a picture (array).
#'
#' @param img [array] a picture loaded by \code{png::readPNG()}.
#' @param alpha [numeric] transparency factor between 0 (transparent) and 1 (opaque).
#' @param color [array] a color name (e.g. "red", "steelblue") or an hexadecimal value (e.g. "#ff0000")
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#' @author Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @export
#'
#' @return
#' This function returns a picture matrix with colors and/or transparency altered.
#' This matrix can be plotted with the function \code{graphics::rasterImage()} (see examples).
#'
#' @examples
#'
#' # Download Rlogo in png format
#' download.file(
#'   url      = "https://www.r-project.org/logo/Rlogo.png",
#'   destfile = "Rlogo.png"
#' )
#'
#' # Import Rlogo (require `png` package)
#' rlogo <- png::readPNG("Rlogo.png")
#'
#' # Get image dimensions
#' xmax <- dim(rlogo)[2]
#' ymax <- dim(rlogo)[1]
#'
#' # Remove figure margins
#' o_par <- par()
#' par(mar = rep(0, 4))
#'
#'
#' ### ORIGINAL IMAGE   ------
#'
#'
#' # Empty plot
#' xmin <- 0
#' ymin <- xmin
#'
#' plot(
#'   x    = xmin,
#'   x    = ymin,
#'   xlim = c(xmin, xmax), # X-axis extent (number of image rows)
#'   ylim = c(ymin, ymax), # Y-axis extent (number of image cols)
#'   asp  = 1,             # Preserve x-axis/y-axis ratio
#'   ann  = FALSE,         # Remove figure annotations
#'   axes = FALSE,         # Remove figure axes
#'   type = "n"            # Remove figure data
#' )
#'
#' # Add Rlogo
#' graphics::rasterImage(
#'   image   = rlogo,
#'   xleft   = xmin,
#'   ybottom = ymin,
#'   xright  = xmax,
#'   ytop    = ymax,
#'   angle   = 0
#' )
#'
#'
#' ### IMAGE WITH TRANSPARENCY   ------
#'
#' # Empty plot
#' plot(
#'   x    = 0,
#'   xlim = c(0, xmax),
#'   ylim = c(0, ymax),
#'   asp  = 1,
#'   ann  = FALSE,
#'   axes = FALSE,
#'   type = "n"
#' )
#'
#' # Add transparency
#' rlogo_alpha <- recolorPhylopic(rlogo, alpha = 0.25, color = NULL)
#'
#' # Add Rlogo
#' graphics::rasterImage(
#'   image   = rlogo_alpha,
#'   xleft   = 0,
#'   ybottom = 0,
#'   xright  = xmax,
#'   ytop    = ymax,
#'   angle   = 0
#' )
#'
#'
#' ### MONOCHROMATIC IMAGE   ------
#'
#' # Empty plot
#' plot(
#'   x    = 0,
#'   xlim = c(0, xmax),
#'   ylim = c(0, ymax),
#'   asp  = 1,
#'   ann  = FALSE,
#'   axes = FALSE,
#'   type = "n"
#' )
#'
#' # Add transparency
#' rlogo_blue <- recolorPhylopic(rlogo, alpha = 1, color = "steelblue")
#'
#' # Add Rlogo
#' graphics::rasterImage(
#'   image   = rlogo_blue,
#'   xleft   = 0,
#'   ybottom = 0,
#'   xright  = xmax,
#'   ytop    = ymax,
#'   angle   = 0
#' )
#'
#' # Delete downloaded Rlogo
#' file.remove("Rlogo.png")
#'
#' # Reset default par()
#' par(o_par)
#'



recolorPhylopic <- function(img, alpha = 0.2, color = NULL) {

  if (!is.null(color)) {

    if (length(color) > 1) {

      stop("`color` must be of length of 1.")

    }

    if (!is.character(color)) {

      stop("`color` must be a color name (e.g. \"red\") or a hexadecimal value.")

    }
  }

  if (!is.null(alpha)) {

    if (!is.numeric(alpha)) {

      stop("`alpha` must be numeric.")

    }

    if (alpha < 0 || alpha > 1) {

      stop("`alpha` must be >= 0 and <= 1.")
    }
  }


  if (is.null(color)) {

    img <- matrix(
      rgb(
        img[ , , 1],
        img[ , , 2],
        img[ , , 3],
        img[ , , 4] * alpha
      ),
      nrow = dim(img)[1]
    )

  } else {

    cols   <- grDevices::col2rgb(color)
    imglen <- length(img[ , , 1])

    img    <- matrix(
      ifelse(
        img[ , , 4] > 0,
        rgb(
          rep(cols[1, 1], imglen),
          rep(cols[2, 1], imglen),
          rep(cols[3, 1], imglen),
          img[ , , 4] * 255 * alpha,
          maxColorValue = 255
        ),
        rgb(
          rep(1, imglen),
          rep(1, imglen),
          rep(1, imglen),
          img[ , , 4] * alpha
        )
      ),
      nrow = dim(img)[1]
    )
  }

  return(img)
}
