#' @title Add inset map to an existing map
#'
#' @description
#' For an internal use only.
#'
#' @param x [raster] a RasterLayer to crop and add into an inset.
#' @param region [vector] the extent of the crop box (xmn, xmx, ymn, ymx).
#' @param where [vector] the coordinates of left-bottom corner of the inset.
#' @param zoom [numeric] a zoom factor.
#' @param title [string] the inset title.
#' @param case [integer] predefined case.
#' @param yat [numeric] a factor to correct title position.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#'
#' @export
#'
#' @return
#' For an internal use only.
#'
#' @examples
#'
#' # No example.
#'



addInset <- function(x = NULL, region = NULL, where = NULL, zoom = 2, title = NULL, case = 1, yat = 0.1){

  region <- as(extent(region), "SpatialPolygons")
  inset  <- crop(x, region)
  box    <- as(extent(inset), "SpatialPolygons")

  xmin <- where[1]
  ymin <- where[2]
  xmax <- xmin + (extent(inset)[2] - extent(inset)[1]) * zoom
  ymax <- ymin + (extent(inset)[4] - extent(inset)[3]) * zoom

  extent(inset) <- c(xmin, xmax, ymin, ymax)

  raster::plotRGB(
    inset, r = 1, g = 2, b = 3,
    colNA   = "#c0e8f3",
    axes    = FALSE,
    add     = TRUE
  )


  if (case == 1) {

    lines(
      x    = rep(xmin + (xmax - xmin) / 2, 2),
      y    = c(extent(box)[3] + (extent(box)[4] - extent(box)[3]) / 2, ymax),
      col  = "white",
      lwd  = 3
    )

    lines(
      x    = c(xmin + (xmax - xmin) / 2, extent(box)[1]),
      y    = rep(extent(box)[3] + (extent(box)[4] - extent(box)[3]) / 2, 2),
      col  = "white",
      lwd  = 3
    )
  }

  if (case == 2) {

    lines(
      x    = rep(xmin + (xmax - xmin) / 10, 2),
      y    = c(extent(box)[3] + (extent(box)[4] - extent(box)[3]) / 2, ymax),
      col  = "white",
      lwd  = 3
    )

    lines(
      x    = c(xmin + (xmax - xmin) / 10, extent(box)[1]),
      y    = rep(extent(box)[3] + (extent(box)[4] - extent(box)[3]) / 2, 2),
      col  = "white",
      lwd  = 3
    )
  }

  if (case == 3) {

    lines(
      x    = rep(extent(box)[1] + (extent(box)[2] - extent(box)[1]) / 2, 2),
      y    = c(extent(box)[3], ymin + (ymax - ymin) / 3),
      col  = "white",
      lwd  = 3
    )

    lines(
      x    = c(xmax, extent(box)[1] + (extent(box)[2] - extent(box)[1]) / 2),
      y    = rep(ymin + (ymax - ymin) / 3, 2),
      col  = "white",
      lwd  = 3
    )
  }

  if (case == 4) {

    lines(
      x    = rep(extent(box)[1] + (extent(box)[2] - extent(box)[1]) / 2, 2),
      y    = c(extent(box)[4], ymin + (ymax - ymin) / 2),
      col  = "white",
      lwd  = 3
    )

    lines(
      x    = c(xmin, extent(box)[1] + (extent(box)[2] - extent(box)[1]) / 2),
      y    = rep(ymin + (ymax - ymin) / 2, 2),
      col  = "white",
      lwd  = 3
    )
  }

  sp::plot(box, add = TRUE, border = "white", lwd = 3, col = "transparent")


  # plot(ocean, col = "#c0e8f3", border = NA, add = TRUE)

  rect(
    xleft    = xmin,
    ybottom  = ymin,
    xright   = xmax,
    ytop     = ymax,
    col      = NA,
    border   = "white",
    lwd      = 3
  )

  if (case == 1) {

    lines(
      x    = rep(xmin + (xmax - xmin) / 2, 2),
      y    = c(extent(box)[3] + (extent(box)[4] - extent(box)[3]) / 2, ymax),
      col  = par()$col,
      lwd  = 1.5
    )

    lines(
      x    = c(xmin + (xmax - xmin) / 2, extent(box)[1]),
      y    = rep(extent(box)[3] + (extent(box)[4] - extent(box)[3]) / 2, 2),
      col  = par()$col,
      lwd  = 1.5
    )
  }

  if (case == 2) {

    lines(
      x    = rep(xmin + (xmax - xmin) / 10, 2),
      y    = c(extent(box)[3] + (extent(box)[4] - extent(box)[3]) / 2, ymax),
      col  = par()$col,
      lwd  = 1.5
    )

    lines(
      x    = c(xmin + (xmax - xmin) / 10, extent(box)[1]),
      y    = rep(extent(box)[3] + (extent(box)[4] - extent(box)[3]) / 2, 2),
      col  = par()$col,
      lwd  = 1.5
    )
  }

  if (case == 3) {

    lines(
      x    = rep(extent(box)[1] + (extent(box)[2] - extent(box)[1]) / 2, 2),
      y    = c(extent(box)[3], ymin + (ymax - ymin) / 3),
      col  = par()$col,
      lwd  = 1.5
    )

    lines(
      x    = c(xmax, extent(box)[1] + (extent(box)[2] - extent(box)[1]) / 2),
      y    = rep(ymin + (ymax - ymin) / 3, 2),
      col  = par()$col,
      lwd  = 1.5
    )
  }

  if (case == 4) {

    lines(
      x    = rep(extent(box)[1] + (extent(box)[2] - extent(box)[1]) / 2, 2),
      y    = c(extent(box)[4], ymin + (ymax - ymin) / 2),
      col  = par()$col,
      lwd  = 1.5
    )

    lines(
      x    = c(xmin, extent(box)[1] + (extent(box)[2] - extent(box)[1]) / 2),
      y    = rep(ymin + (ymax - ymin) / 2, 2),
      col  = par()$col,
      lwd  = 1.5
    )
  }

  sp::plot(box, add = TRUE, border = par()$col, lwd = 1.5, col = "transparent")

  rect(
    xleft    = xmin,
    ybottom  = ymin,
    xright   = xmax,
    ytop     = ymax,
    col      = NA,
    border   = par()$col,
    lwd      = 1.5
  )

  text(
    x       = xmax,
    y       = ymax - yat * (ymax - ymin),
    labels  = title,
    col     = par()$col.axis,
    font    = 2,
    pos     = 2
  )
}
