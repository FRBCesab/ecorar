#' Panel of Two Scatterplots of PCoA Axis
#'
#' This script produces the Loiseau, Mouquet et al.'s 2020 paper Figure 2, i.e.
#' a panel of two scatterplots (mammals on the left and birds on the right)
#' representing functional spaces based on species functional traits.
#'
#' @author Nicolas Casajus, \email{nicolas.casajus@@fondationbiodiversite.fr},
#'         Nicolas Loiseau, \email{nicolas.loiseau1@@gmail.com},
#'         Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr}
#'
#' @date 2020/06/10


## Parameters ----

n     <- 1                          # ID of the first sub-plot
plots <- list()                     # Sub-plots storage


for (taxa in taxas) {


  ## Subset Data for Taxa ----

  subdatas <- datas[datas$"class" == taxa, ]
  subdatas <- subdatas[ , c(paste0("pcoa_axis_", pcoa_axes[[taxa]]), "din", "dr_class")]

  subdatas <- subdatas[!is.na(subdatas$"dr_class"), ]

  colnames(subdatas)[1:2] <- c("x", "y")
  subdatas$"x" <- jitter(subdatas$"x", jitter_val)
  subdatas$"y" <- jitter(subdatas$"y", jitter_val)


  ## Add Scatterplot ----

  gplot <- ggplot(subdatas, aes(x, y, fill = din)) +
    
    geom_point(shape = 21, size = 3, colour = "transparent") +
    
    scale_fill_gradientn(colours = color_distinctiveness)

  
  ## Get Plot Area Coordinates ----
  
  x_lims <- ggplot_build(gplot)$layout$panel_params[[1]]$x.range
  lag_x  <- diff(x_lims) * 0.04
  y_lims <- ggplot_build(gplot)$layout$panel_params[[1]]$y.range
  lag_y  <- diff(y_lims) * 0.04
  
  
  ## Add Convex Hulls ----

  for (classe in classes) {

    coords <- subdatas[subdatas$"dr_class" == classe, ]
    hull   <- findHull(coords[complete.cases(coords), ])

    gplot <- gplot +
      
      geom_polygon(data = hull, alpha = 0.1, size = 2.0, fill = "white",
                   colour = "white") +
      
      geom_polygon(data = hull, alpha = 0.1, size = 1.0, fill = "transparent",
                   colour = color_classes[classe])
  }


  ## Change Theme ----

  gplot <- gplot +
    
    labs(x = paste("PCoA axis", pcoa_axes[[taxa]][1]),
         y = paste("PCoA axis", pcoa_axes[[taxa]][2])) +
    
    gg_theme


  ## Add Din Legend (right panel) ----

  if (taxa == "mammals") {

    gplot <- gplot + 
      
      theme(legend.position = "None")
  }

  if (taxa == "birds") {

    gplot <- gplot +

      theme(legend.position = c(0.88, 0.88)) +
      
      labs(fill = "Distinctiveness")
  }


  ## Add DR Class Legend (left panel) ----

  if (taxa == "mammals") {

    segs <- data.frame(
      x     = rep(0.21, 3),
      xend  = rep(0.25, 3),
      y     = seq(0.47, 0.42, by = -0.025),
      label = classes_labs
    )

    gplot <- gplot +

      geom_rect(aes(xmin = segs[1, "x"] - 0.02, 
                    xmax = segs[1, "xend"] + 0.16, 
                    ymin = 0.40, 
                    ymax = 0.49),
                fill = "white", color = "white", size = 0.25)

    for (i in 1:nrow(segs)) {

      gplot <- gplot +

      annotate(geom = "segment", x = segs[i, "x"], y = segs[i, "y"],
               xend = segs[i, "xend"], yend = segs[i, "y"], size = 1,
               color = color_classes[i]) +

      annotate(geom = "text", x = segs[i, "xend"] + 0.01, y = segs[i, "y"], 
               label = segs[i, "label"], hjust = "left", size = 8.5,
               family  = par_family, colour = dark_grey, fontface = 1)
    }
  }


  ## Add Sub-Plot Label ----
  
  label <- ifelse(taxa == "mammals", "a", "b")
    
  gplot <- gplot +

    annotate(geom = "text", x = x_lims[1] + lag_x, y = y_lims[2] - lag_y, 
             label = label, size = 12, family  = par_family, fontface = 2, 
             colour = dark_grey)


  ## Add Taxa Silhouette ----

  gplot <- cowplot::ggdraw(gplot) +

    cowplot::draw_image(image = icons[[taxa]],
                        x     = ifelse(taxa == "birds", -0.32, -0.25),
                        y     = ifelse(taxa == "birds",  0.40,  0.40),
                        scale = ifelse(taxa == "birds",  0.12,  0.15))

  plots[[n]] <- gplot
  n          <- n + 1

}


## Arrange Sub-Plots ----

mat <- matrix(c(1, 2, 1, 2), ncol = 2, nrow = 2, byrow = TRUE)

grobs <- gridExtra::arrangeGrob(plots[[1]], plots[[2]], layout_matrix = mat)


## Export Figure ----

ggsave(
  filename  = here::here("figures", paste0(figname, ".png")),
  plot      = grobs,
  width     = 24,
  height    = 12,
  units     = "in",
  dpi       = 600
)


## Message ----

usethis::ui_done(
  paste(
    usethis::ui_value(figname), 
    "successfully exported in figures/"
  )
)
