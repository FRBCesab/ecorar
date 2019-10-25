#' --------------------------------------------------------------------------   @Header
#'
#' @title Figure 1 - Panel of Pcoa scatterplots
#'
#' @description
#' This R script...
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @date 2019/08/05
#'
#' --------------------------------------------------------------------------   @VVVVVV




#'  -------------------------------------------------------------------------   @Parameters


n        <- 1                          # ID of first plot
plots    <- list()                     # Subplots storage



for (taxa in taxas) {



#' ---------------------------------------------------------------------------- @SubsetData


  subdatas <- datas[datas[ , "class"] == taxa, ]
  subdatas <- subdatas[ , c(paste0("pcoa_axis_", pcoa_axes[[taxa]]), "din", "dr_class")]

  subdatas <- subdatas[!is.na(subdatas[ , "dr_class"]), ]
  # subdatas <- subdatas[subdatas[ , "dr_class"] %in% classes, ]

  colnames(subdatas)[1:2] <- c("x", "y")
  subdatas[ , "x"] <- jitter(subdatas[ , "x"], jitter_val)
  subdatas[ , "y"] <- jitter(subdatas[ , "y"], jitter_val)



#' ---------------------------------------------------------------------------- @PlotSpecies


  gplot <- ggplot(
    data     = subdatas,
    mapping  = aes(x, y, fill = din)
  ) +

  geom_point(
    shape   = 21,
    size    = 3,
    colour  = "transparent"
  ) +

  scale_fill_gradientn(colours = color_distinctiveness)



#' ---------------------------------------------------------------------------- @AddHulls


  for (classe in classes) {

    coords <- subdatas[subdatas[ , "dr_class"] == classe, ]
    coords <- coords[complete.cases(coords), ]
    hull   <- findHull(coords)

    gplot <- gplot +

    geom_polygon(
      data    = hull,
      alpha   = 0.1,
      size    = 2.0,
      colour  = "white",
      fill    = "white"
    ) +

    geom_polygon(
      data    = hull,
      alpha   = 0.1,
      size    = 1.0,
      colour  = color_classes[classe],
      fill    = "transparent"
    )
  }



#' ---------------------------------------------------------------------------- @PlotTheme


  gplot <- gplot +

  labs(
    x = paste("PCoA axis", pcoa_axes[[taxa]][1]),
    y = paste("PCoA axis", pcoa_axes[[taxa]][2])
  ) +

  theme_bw() +

  theme(
    text        = element_text("serif"),
    axis.text   = element_text(size = 18),
    axis.title  = element_text(size = 20)
  )



#' ---------------------------------------------------------------------------- @AddDinLegend


  if (taxa == "mammals") {

    gplot <- gplot + theme(legend.position = "None")
  }

  if (taxa == "birds") {

    gplot <- gplot +

    theme(
      legend.position    = c(0.88, 0.88),
      legend.title       = element_text(size = 16),
      legend.text        = element_text(size = 14),
      legend.background  = element_rect(fill = "white", color = light_grey)
    ) +

    labs(fill = "Distinctiveness")
  }



#' ---------------------------------------------------------------------------- @AddDRClassLegend


  if (taxa == "mammals") {

    segs <- data.frame(
      x     = rep(0.12, 3),
      xend  = rep(0.16, 3),
      y     = seq(0.47, 0.43, by = -0.02),
      label = classes_labs
    )

    gplot <- gplot +

    geom_rect(
      mapping  = aes(
        xmin = 0.11, xmax = 0.25, ymin = 0.42, ymax = 0.48
      ),
      fill     = "white",
      color    = light_grey,
      size     = 0.25
    )

    for (i in 1:3) {

      gplot <- gplot +

      annotate(
        geom    = "segment",
        x       = segs[i, "x"],
        y       = segs[i, "y"],
        xend    = segs[i, "xend"],
        yend    = segs[i, "y"],
        color   = color_classes[i],
        size    = 1
      ) +

      annotate(
        geom    = "text",
        x       = segs[i, "xend"] + 0.01,
        y       = segs[i, "y"],
        label   = segs[i, "label"],
        hjust   = "left",
        size    = 6,
        family  = "serif"
      )

    }
  }



#' ---------------------------------------------------------------------------- @AddSilhouette


  gplot <- cowplot::ggdraw(gplot) +

  cowplot::draw_image(
    image  = icons[[taxa]],
    x      = ifelse(taxa == "birds", -0.32, -0.30),
    y      = ifelse(taxa == "birds",  0.40,  0.40),
    scale  = ifelse(taxa == "birds",  0.12,  0.15)
  )



  plots[[n]] <- gplot

  n <- n + 1

}



#' ---------------------------------------------------------------------------- @ArrangeSubplots


mat <- matrix(
  data   = c(1, 2, 1, 2),
  ncol   = 2,
  nrow   = 2,
  byrow  = TRUE
)

grobs <- gridExtra::arrangeGrob(
  plots[[1]], plots[[2]],
  layout_matrix = mat
)



#' ---------------------------------------------------------------------------- @ExportFigure


ggsave(
  filename  = file.path(path_figs, paste0(figname1, ".png")),
  plot      = grobs,
  width     = 24,
  height    = 12,
  units     = "in",
  dpi       = 600
)
