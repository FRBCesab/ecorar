compassRose <- function (x, y = x, rot = 0, cex.cr = 1, cex.let = 1,
  labels = c("S", "W", "N", "E"), offset = 1.2, col = c("#95D8EB", "#888888"),
  border = c("#888888", "#888888"), lty = 1, col.let = "#888888") {

  wh  <- graphics::strheight("M")
  rot <- pi * rot/180

  mat.rot <- matrix(c(cos(rot), sin(rot), -sin(rot), cos(rot)), 2)

  lwh <- cex.cr * 4.00 * wh
  swh <- cex.cr * 0.45 * wh

  rex <- rep(c(0, -1, 0, 1), each = 2)
  rey <- rep(c(-1, 0, 1, 0), each = 2)

  rex1 <- c(lwh * rex, swh * c(-1, 1) * rey + rex * swh)
  rey1 <- c(lwh * rey, swh * c(-1, 1) * rex + rey * swh) * rep(c(1, -1), each = 2, 4)

  matxy <- as.matrix(cbind(rex1, rey1))
  matxy <- matxy %*% mat.rot

  cr.col <- rep(col, length.out = 8)
  cr.bd  <- rep(border, length.out = 8)

  for (i in 1:8) {
      graphics::polygon(
        x + c(0, matxy[i, 1], matxy[8 + i, 1]),
        y + c(0, matxy[i, 2], matxy[8 + i, 2]),
        col = cr.col[i], border = cr.bd[i], lty = lty
      )
  }

  graphics::text(
    x + offset * matxy[seq(1, by = 2, length.out = 4), 1],
    y + offset * matxy[seq(1, by = 2, length.out = 4), 2],
    labels, cex = cex.let, col = col.let)


  invisible(NULL)
}
