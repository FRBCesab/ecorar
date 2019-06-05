path_data <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/data/"           ### !!!
path_R    <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/RALLL/scripts/"  ### !!!

source(paste0(path_R, "__recolorPhylopic.R"))
source(paste0(path_R, "__addPhylopic.R"))


cols <- c("#00AFBB", "#E7B800", "orangered")


load(paste0(path_data, "threats.RData"))

datas <- threats

datas$Threats  <- as.character(datas$Threats)
datas$Taxa  <- as.character(datas$Taxa)
datas$DR_class <- as.character(datas$DR_class)


values_m  <- tapply(datas$Value, list(datas$Threats, datas$Taxa, datas$DR_class), function(x) mean(x, na.rm = TRUE))
values_sd <- tapply(datas$Value, list(datas$Threats, datas$Taxa, datas$DR_class), function(x) sd(x, na.rm = TRUE))


iucn <- list()
iucn[[1]] <- as.data.frame(matrix(c(1853, 483, 331, 410, 69, 22, 132, 3, 220), ncol = 3, byrow = TRUE))
rownames(iucn[[1]]) <- c("LC", "NE", "TH")
colnames(iucn[[1]]) <- c("AVG", "D25R25", "D75R75")
names(iucn)[1] <- "birds"

iucn[[2]] <- as.data.frame(matrix(c(907, 196, 32, 103, 0, 39, 116, 4, 166), ncol = 3, byrow = TRUE))
rownames(iucn[[2]]) <- c("LC", "NE", "TH")
colnames(iucn[[2]]) <- c("AVG", "D25R25", "D75R75")
names(iucn)[2] <- "mammals"

for (i in 1:2)
  iucn[[i]] <- apply(iucn[[i]], 2, function(x) x / sum(x) * 100)

labels <- c("Climate change", "Mean conflict", "Mean HDI", "Target achievement (%)")

icons <- list()
icons[["mammals"]] <- image_data("5a5dafa2-6388-43b8-a15a-4fd21cd17594", size = 512)[[1]]
icons[["mammals"]] <- icons[["mammals"]][ , dim(icons[["mammals"]])[2]:1, ]
icons[["birds"]]   <- image_data("c3c19b65-cb8d-40bd-b1a6-82a3277bcd4f", size = 512)[[1]]



png(
  file       = paste0("~/Desktop/figure-3.png"),
  width      = 12.00,
  height     =  7.00,
  units      = "in",
  res        = 600,
  pointsize  = 18
)

par(
  # xaxs      = "i",
  # yaxs      = "i",
  family    = "serif",
  mar       = rep(1, 4),
  cex.axis  = 1.25,
  mgp       = c(2, .5, 0),
  tcl       = -0.25,
  xpd       = FALSE,
  new       = FALSE,
  fig       = c(0, 1, 0, 1),
  col       = "#666666",
  col.axis  = "#666666",
  cex.axis  = 0.75,
  fg        = "#666666"
)
par(mfcol = c(2, 5))


comp <- 1


for (j in 1:dim(values_m)[1]) {

  for (i in 1:dim(values_m)[2]) {

    # if (comp / 2 == round(comp / 2)) {
    #   if (comp == 2) {
    #     par(mar = c(1, 1, 0, 0))
    #   } else {
    #     par(mar = c(.5, 1, .5, 0))
    #   }
    # } else {
    #   if (comp == 1) {
    #     par(mar = c(.5, 1, .5, 0))
    #   } else {
    #     par(mar = c(.5, .5, .5, .5))
    #   }
    # }

    if (comp / 2 == round(comp / 2)) {
      par(mar = c(1, 1, 0.5, 0.5))
    } else {
      par(mar = c(0.5, 1, 1, .5))
    }

    plot(
      0,
      xlim = c(
        min(values_m[j, , ] - values_sd[j, , ], na.rm = TRUE),
        max(values_m[j, , ] + values_sd[j, , ], na.rm = TRUE)
      ),
      ylim = c(0.5, 3.5),
      axes = FALSE,
      type = "n",
      ann = FALSE
    )

    grid()

    if (!is.na(values_m[j, i, 1])) {

      for (k in 1:3){
        lines(
          x = c(
            values_m[j, i, k] - values_sd[j, i, k],
            values_m[j, i, k] + values_sd[j, i, k]
          ),
          y = rep(k, 2),
          col = cols[k],
          lwd = 2
        )
        lines(
          x = rep(
            values_m[j, i, k] - values_sd[j, i, k],
            2
          ),
          y = c(k-.1, k + .1),
          col = cols[k],
          lwd = 2
        )
        lines(
          x = rep(
            values_m[j, i, k] + values_sd[j, i, k],
            2
          ),
          y = c(k-.1, k + .1),
          col = cols[k],
          lwd = 2
        )
        points(x = values_m[j, i, k], y = k, col = cols[k], pch = 19, cex = 2)
      }
    }
    box()

    if (comp < 3) {
      par(mgp = c(.5, 0, 0))
      axis(2, at = 1:3, labels = dimnames(values_m)[3][[1]], lwd = 0, font = 2)
    }

    if (comp / 2 == round(comp / 2)) {
      par(mgp = c(.5, -.05, 0))
      axis(1, lwd = 0, font = 2)
    }

    if (comp / 2 != round(comp / 2)) {
      mtext(text = labels[j], side = 3, font = 2, cex = .85)
    }

    if (comp == 1) {

      addPhylopic(
        img    = icons[[2]],
        x      = .80,
        y      = .90,
        ysize  = .20,
        alpha  = 1,
        color  = "#777777",
        AR     = .5
      )
    }

    if (comp == 2) {

      addPhylopic(
        img    = icons[[1]],
        x      = .80,
        y      = .90,
        ysize  = .20,
        alpha  = 1,
        color  = "#777777",
        AR     = .5
      )
    }

    comp <- comp + 1
  }
}

cols <- c("#026666", "#aaaaaa", "#c53131")
sloc <- c("#f7f7f7", "#333333", "#F6c8c8")

for (i in 1:2) {

  if (i == 2) {
    par(mar = c(1, 1, 0.5, 0.5))
  } else {
    par(mar = c(0.5, 1, 1, .5))
  }

  plot(
    0,
    xlim = c(
      -10,
      110
    ),
    ylim = c(0.5, 3.5),
    axes = FALSE,
    type = "n",
    ann = FALSE
  )

  grid()
  box()

  for (j in 1:3) {

     # pos <- order(iucn[[i]][ , j], decreasing = TRUE)
     # for (k in pos) {
     #   lines(x = c(0, iucn[[i]][k, j]), y = rep(j, 2), col = cols[k], lwd = 2)
     # }

    for (k in c(2, 1, 3)) {
      points(iucn[[i]][k, j], j, cex = 4, col = cols[k], pch = 19)
      text(iucn[[i]][k, j], j, labels = rownames(iucn[[i]])[k], cex = .85, col = sloc[k], family = "sans", font = 2)
    }
  }

  if (i == 2) {
    par(mgp = c(.5, -.05, 0))
    axis(1, lwd = 0, font = 2)
  } else {
    mtext(text = "IUCN status (%)", side = 3, font = 2, cex = .85)
  }
}

dev.off()
