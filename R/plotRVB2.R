#' Apply a color gradient to a raster
#'
#' This function applies a color gradient to a raster where color is defined by the user.
#'
#' @param x [RasterLayer] a RasterLayer to plot.
#' @param n_classes [...] ...
#' @param breaks [...] ...
#' @param palettes [...] ...
#' @param alpha [integer] transparency factor between 0 (transparent) and 255 (opaque).
#' @param bgalpha [integer] transparency factor for NA values between 0 (transparent) and 255 (opaque).
#' @param colNA [string] color for the background (i.e. NA values).
#' @param axes [boolean] If TRUE, axes are added.
#' @param xlab [string] x-axis label.
#' @param ylab [string] y-axis label.
#' @param add [boolean] If TRUE, raster is added the active graphical device.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#'
#' @export
#'
#' @return
#' This function returns a 3-bands raster with the first layer is for the Red,
#' the second for the Green and the third for the Blue. Each non NA cell contains
#' the quantity of a specific color (R, G or B).
#'
#' @examples
#'
#' # No example.



plotRVB2 <- function(x, n_classes = NULL, breaks = NULL, palettes = NULL,
  alpha = 255, bgalpha = 0, colNA = NA, axes = FALSE, xlab = "", ylab = "",
  add = FALSE) {


  ## Argument checks ----

  if (missing(x)) {
    stop("Missing `x`: single spatial raster layer is required.")
  }

  if (is.null(x)) {
    stop("Single spatial raster layer is required.")
  }

  if (class(x) != "RasterLayer") {
    stop("Single spatial raster layer is required.")
  }


  ## Identify non-NA cells ----

  cells <- which(!is.na(x[]))

  classes <- list(
    which(x[][cells] <= breaks[1]),
    which(x[][cells] >= breaks[1] & x[][cells] <= breaks[2]),
    which(x[][cells] >= breaks[2])
  )

  diff_1 <- diff(range(x[][cells[classes[[1]]]]))
  diff_2 <- diff(range(x[][cells[classes[[3]]]]))
  
  diff_1 <- ifelse(is.finite(diff_1), diff_1, NA)
  diff_2 <- ifelse(is.finite(diff_2), diff_2, NA)
  
  class_by <- sum(c(diff_1, diff_2), na.rm = TRUE) / n_classes

  if (!is.na(diff_1)) {
    cols_1 <- colorRampPalette(palettes[[1]])(
      length(
        seq(
          min(range(x[][cells[classes[[1]]]])),
          breaks[1],
          by = class_by
        )
      ) - 1
    )
  } else {
    cols_1 <- NULL
  }
  
  if (!is.na(diff_2)) {
    cols_2 <- colorRampPalette(palettes[[3]])(
      length(
        seq(
          breaks[2],
          max(range(x[][cells[classes[[3]]]])),
          by = class_by
        )
      ) - 1
    )
  } else {
    cols_2 <- NULL
  }
  
  colors <- list(cols_1, palettes[[2]], cols_2)

  red <- green <- blue <- x

  for (i in 1:length(colors)) {

    if (length(classes[[i]])) {
      
      x1 <- x
      x1[][-cells[classes[[i]]]] <- NA
      
      hexa <- leaflet::colorNumeric(
        palette  = colors[[i]],
        domain   = values(x1),
        na.color = NA
      )
      hexa <- hexa(values(x1))
      rgb  <- grDevices::col2rgb(hexa)
      
      raster::values(red)[cells[classes[[i]]]]   <- rgb[1, cells[classes[[i]]]]
      raster::values(green)[cells[classes[[i]]]] <- rgb[2, cells[classes[[i]]]]
      raster::values(blue)[cells[classes[[i]]]]  <- rgb[3, cells[classes[[i]]]]
    }
  }

  ras <- stack(red, green, blue)


  ## Plot the raster ----

  plotRGB(
    ras, r = 1, g = 2, b = 3,
    alpha   = alpha,
    bgalpha = bgalpha,
    colNA   = colNA,
    axes    = axes,
    xlab    = xlab,
    ylab    = ylab,
    add     = add
  )


  ## Add Legend Colors ----

  ybottom <- -6550000
  xstart  <- -16500000
  xleft   <- xstart

  ccolors <- colors
  colors <- unlist(colors)

  for (i in 1:length(colors)) {

    rect(
      xleft    = xleft,
      ybottom  = ybottom,
      xright   = xleft + 35000,
      ytop     = ybottom + 350000,
      border   = NA,
      col      = colors[i]
    )

    xleft   <- xleft + 35000
  }


  ## Add Legend Box ----

  rect(
    xleft    = xstart,
    ybottom  = ybottom,
    xright   = xleft,
    ytop     = ybottom + 350000,
    border   = "white",
    col      = NA
  )



#'  -------------------------------------------------------------------------   @AddLegendLabel


  text(
    x       = xstart,
    y       = ybottom + 200000,
    labels  = round(min(x[][cells]), 2),
    pos     = 1,
    col     = par()$col.axis,
    font    = 2,
    cex     = .65
  )

  text(
    x       = xleft,
    y       = ybottom + 200000,
    labels  = round(max(x[][cells]), 2),
    pos     = 1,
    col     = par()$col.axis,
    font    = 2,
    cex     = .65
  )

  if (!is.na(diff_1)) {
    
    text(
      x       = xstart + (35000 * (length(ccolors[[1]]) + 1)),
      y       = ybottom + 200000,
      labels  = 0,
      pos     = 1,
      col     = par()$col.axis,
      font    = 2,
      cex     = .65
    )
  }

  text(
    x       = xstart + ((xleft - xstart) / 2),
    y       = ybottom + 200000,
    labels  = "Null Model",
    pos     = 3,
    col     = par()$col.axis,
    font    = 2,
    cex     = .85
  )


  return(ras)
}
