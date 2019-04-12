addGraticules <- function(
  prj           = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs",
  parallels     = seq(-60,  60, by =  30),
  meridians     = seq(-160, 160, by = 40),
  line.color    = "#aaaaaa",
  line.size     = 0.50,
  line.type     = 1,
  labels        = TRUE,
  lang          = "en",
  text.color    = "#555555",
  text.size     = 0.45,
  text.pos      = c(0, 0, 1),
  text.rotation = TRUE,
  add           = FALSE

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

      parallels[[i]] <- Lines(
        list(
          Line(
            cbind(
              xext,
              rep(lats[i], length(xext))
            )
          )
        ),
        as.character(i)
      )
    }

    parallels <- SpatialLinesDataFrame(
      sl   = SpatialLines(
        LinesList   = parallels,
        proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      ),
      data = data.frame(latitude = lats)
    )
    parallels <- spTransform(parallels, CRSobj = prj)


    ### PLOT PARALLELS --------------------------

    if (add) {
      plot(
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

      meridians[[i]] <- Lines(
        list(
          Line(
            cbind(
              rep(lons[i], length(yext)),
              yext
            )
          )
        ),
        as.character(i)
      )
    }

    meridians <- SpatialLinesDataFrame(
      sl   = SpatialLines(
        LinesList   = meridians,
        proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      ),
      data = data.frame(longitude = lons)
    )
    meridians <- spTransform(meridians, CRSobj = prj)


    ### PLOT MERIDIANS --------------------------

    if (add) {
      plot(
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
