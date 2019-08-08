addPhylopic <- function(img, x = NULL, y = NULL, ysize = NULL, alpha = 0.2, color = NULL, AR = 1) {

  img  <- recolorPhylopic(img, alpha, color)
  dims <- dim(img)[1:2]

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
