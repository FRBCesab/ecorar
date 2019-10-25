#' @title Apply a color gradient to a raster
#'
#' @description
#' This function applies a color gradient to a raster where color is defined by the user.
#'
#' @param x [RasterLayer] a RasterLayer to plot.
#' @param zlim [numeric] a vector of length 2. Range of values to plot.
#' @param type [string] not implemented.
#' @param palette [string] not well implemented.
#' @param reverse [boolean] not well implemented.
#' @param alpha [integer] transparency factor between 0 (transparent) and 255 (opaque).
#' @param bgalpha [integer] transparency factor for NA values between 0 (transparent) and 255 (opaque).
#' @param colNA [string] color for the background (i.e. NA values).
#' @param axes [boolean] If TRUE, axes are added.
#' @param xlab [string] x-axis label.
#' @param ylab [string] y-axis label.
#' @param add [boolean] If TRUE, raster is added the active graphical device.
#' @param zmax [numeric] not implemented.
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



plotRVB <- function(
  x = NULL, zlim = NULL, type = NULL, palette = "Spectral", reverse = TRUE,
  alpha = 255, bgalpha = 0, colNA = NA, axes = FALSE, xlab = "", ylab = "",
  add = FALSE, zmax = NULL
) {


  if (is.null(x)) {

    stop("Single spatial raster layer is required.")

  }

  if (class(x) != "RasterLayer") {

    stop("Single spatial raster layer is required.")

  }


  ### List of Brewer Palettes                                                   ----------

  rampcolors <- data.frame(
    palette          = rownames(brewer.pal.info),
    maxcolors        = brewer.pal.info[ , "maxcolors"],
    stringsAsFactors = FALSE
  )


  ### Check Palette                                                             ----------

  pos <- which(tolower(rampcolors[ , "palette"]) == tolower(palette))

  if (length(pos) == 0) {

    stop("Wrong palette name. See brewer.pal.info for a list of available palettes.")

  }


  ### Identify non-NA cells                                                     ----------

  cells <- which(!is.na(x[]))


  ### Get colors defining palette (from 8 to 12)                                ----------

  pal <- RColorBrewer::brewer.pal(
    n    = rampcolors[pos, "maxcolors"],
    name = gsub("-rev", "", palette)
  )


  ### Reverse palette (if required)                                             ----------

  if (reverse) pal <- pal[length(pal):1]


  ### Special case: Binary data                                                 ----------

  # if (!is.null(type)) {
  #
  #   PAL <- pal[c(3, length(pal)-2)]
  #
  #   unique_value <- unique(values(x))
  #
  #   if (length(unique_value) == 3) { # Presence & Absence data (+ NA)
  #
  #     pal <- PAL
  #
  #   } else {
  #
  #     if (0 %in% unique_value) { # Only Absence data (+ NA)
  #
  #       pal <- PAL[1]
  #
  #     } else {
  #
  #       pal <- PAL[length(PAL)] # Only Presence data (+ NA)
  #
  #     }
  #   }
  # }


  ### Create colors gradient based on raster values                             ----------

  if (min(unique(values(x)), na.rm = TRUE) == 0) {

    colors <- c(
      "#aaaaaa",
      colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(255)
    )

  } else {

    colors <- colorRampPalette(brewer.pal(name = "YlOrRd", n = 9))(255)

  }

  # if (is.null(zmax)) {

    hexa <- leaflet::colorNumeric(
      palette  = colors,
      domain   = values(x),
      na.color = NA
    )

  # } else {
  #
  #   hexa <- leaflet::colorNumeric(
  #     palette  = colors,
  #     domain   = seq(0, zmax, by = 1),
  #     na.color = NA
  #   )
  # }

  hexas <- hexa(values(x))


  ### Convert Hexa to RBG                                                       ----------

  rgbs  <- grDevices::col2rgb(hexas)


  ### Create a 3-bands RGB RasterStack                                          ----------

  red <- green <- blue <- x

  raster::values(red)[cells]   <- rgbs[1, cells]
  raster::values(green)[cells] <- rgbs[2, cells]
  raster::values(blue)[cells]  <- rgbs[3, cells]

  x <- stack(red, green, blue)


  ### Plot the raster                                                           ----------

  plotRGB(
    x, r = 1, g = 2, b = 3,
    zlim    = zlim,
    alpha   = alpha,
    bgalpha = bgalpha,
    colNA   = colNA,
    axes    = axes,
    xlab    = xlab,
    ylab    = ylab,
    add     = add
  )

  return(x)
}
