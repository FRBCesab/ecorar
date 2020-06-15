#' @title Add graticules (parallels and meridians) to a map
#'
#' @description
#' This function adds graticules (parallels and meridians) to a map and returns
#' axes labels and coordinates.
#'
#' @param prj [string] the CRS of the map in the PROJ4 standard.
#' @param parallels [numeric] vector of latitudes to add parallels.
#' @param meridians [numeric] vector of longitudes to add meridians.
#' @param line.color [string] color of graticules.
#' @param line.size [numeric] width of graticules.
#' @param line.type [integer] type of graticules.
#' @param add [boolean] If TRUE, graticules are added to the map.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#'
#' @export
#'
#' @return
#' This function returns a 2-elements list.
#'
#' A first data frame for the meridians with 4 columns:
#'   - x: the x-coordinates of meridians
#'   - y: the y-coordinates of meridians
#'   - label: the labels of meridians (in degrees)
#'   - direction: the direction, i.e. W(est) or E(ast)
#'
#' A second data frame for the parallels with 4 columns:
#'   - x: the x-coordinates of parallels
#'   - y: the y-coordinates of parallels
#'   - label: the labels of parallels (in degrees)
#'   - direction: the direction, i.e. S(outh) or N(orth)
#'
#' @examples
#'
#' library(maps)
#'
#' maps::map()
#' addGraticules(
#'   prj = "+init=epsg:4326",
#'   add = TRUE
#' )
#'
#' ## [[1]]
#' ##      x      y label direction
#' ## 1 -160 -90.00   160         W
#' ## 2 -120 -89.95   120         W
#' ## 3  -80 -89.90    80         W
#' ## 4  -40 -89.85    40         W
#' ## 5    0 -89.80     0
#' ## 6   40 -89.75    40         E
#' ## 7   80 -89.70    80         E
#' ## 8  120 -89.65   120         E
#' ## 9  160 -89.60   160         E
#'
#' ## [[2]]
#' ##        x   y label direction
#' ## 1 -180.0 -60    60         S
#' ## 2 -179.9 -30    30         S
#' ## 3 -179.8   0     0
#' ## 4 -179.7  30    30         N
#' ## 5 -179.6  60    60         N
#'



addGraticules <- function(
  prj        = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs",
  parallels  = seq( -60,  60, by = 30),
  meridians  = seq(-160, 160, by = 40),
  line.color = "#aaaaaa",
  line.size  = 0.50,
  line.type  = 1,
  add        = FALSE
) {


  if (!is.null(parallels)) {

    if (is.null(dev.list())) {

      stop("No graphical device found.")
    }


    ### CREATE PARALLELS ------------------------

    lats      <- parallels
    xext      <- seq(-180, 180, by = 0.1)
    parallels <- list()

    for (i in 1:length(lats)) {

      parallels[[i]] <- sp::Lines(
        list(
          sp::Line(
            cbind(
              xext,
              rep(lats[i], length(xext))
            )
          )
        ),
        as.character(i)
      )
    }

    parallels <- sp::SpatialLinesDataFrame(
      sl   = sp::SpatialLines(
        LinesList   = parallels,
        proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      ),
      data = data.frame(latitude = lats)
    )
    parallels <- sp::spTransform(parallels, CRSobj = prj)


    ### PLOT PARALLELS --------------------------

    if (add) {
      sp::plot(
        parallels,
        add = TRUE,
        col = line.color,
        lwd = line.size,
        lty = line.type
      )
    }
  }


  if (!is.null(meridians)) {

    if (is.null(dev.list())) {

      stop("No graphical device found.")
    }


    ### CREATE MERIDIANS ------------------------

    lons      <- meridians
    yext      <- seq(-90, 90, by = 0.05)
    meridians <- list()

    for (i in 1:length(lons)) {

      meridians[[i]] <- sp::Lines(
        list(
          sp::Line(
            cbind(
              rep(lons[i], length(yext)),
              yext
            )
          )
        ),
        as.character(i)
      )
    }

    meridians <- sp::SpatialLinesDataFrame(
      sl   = sp::SpatialLines(
        LinesList   = meridians,
        proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      ),
      data = data.frame(longitude = lons)
    )
    meridians <- sp::spTransform(meridians, CRSobj = prj)


    ### PLOT MERIDIANS --------------------------

    if (add) {
      sp::plot(
        meridians,
        add = TRUE,
        col = line.color,
        lwd = line.size,
        lty = line.type
      )
    }
  }


  ### GET PARALLELS LABELS COORDS --------------------

  lat_coords <- data.frame()

  if (!is.null(lats)) {

    for (i in lats) {

      pos <- which(parallels@data$latitude == i)
      xy  <- parallels@lines[[pos]]@Lines[[1]]@coords

      if (i == 0) {
        lab <- ""
      } else {
        if (i > 0) {
          lab <- "N"
        } else {
          lab <- "S"
          i <- -1 * i
        }
      }

      sop <- which.min(xy[ , 2])
      tmp <- data.frame(
        x         = xy[pos, 1],
        y         = xy[pos, 2],
        label     = i,
        direction = lab
      )
      lat_coords <- rbind(lat_coords, tmp)
    }
  }

  lon_coords <- data.frame()

  if (!is.null(lons)) {

    for (i in lons) {

      pos <- which(meridians@data$longitude == i)
      xy  <- meridians@lines[[pos]]@Lines[[1]]@coords

      if (i == 0) {
        lab <- ""
      } else {
        if (i > 0) {
          lab <- "E"
        } else {
          lab <- "W"
          i <- -1 * i
        }
      }

      sop <- which.min(xy[ , 1])
      tmp <- data.frame(
        x         = xy[pos, 1],
        y         = xy[pos, 2],
        label     = i,
        direction = lab
      )
      lon_coords <- rbind(lon_coords, tmp)
    }
  }

  return(
    list(
      lon_coords,
      lat_coords
    )
  )

}
