#' Panel of Two Scatterplots of Distinctiveness and Restrictiveness
#'
#' This script produces the Loiseau, Mouquet et al.'s 2020 paper Figure 1, i.e.
#' a panel of two scatterplots (mammals on the left and birds on the right)
#' representing functional distinctiveness and geographical restrictiveness.
#'
#' @author Nicolas Casajus, \email{nicolas.casajus@@fondationbiodiversite.fr},
#'         Nicolas Loiseau, \email{nicolas.loiseau1@@gmail.com},
#'         Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr}
#'
#' @date 2020/06/10


owarn <- options()$warn
options(warn = -1)


## Parameters ----

n     <- 1                          # ID of the first sub-plot
plots <- list()                     # Sub-plots storage



for (taxa in taxas) {
  
  
  ## Subset Data for Taxa ----
  
  subdatas <- datas[datas$"class" == taxa, ]
  subdatas <- subdatas[ , c( "din", "rin","dr_class")]

  
  ## Get Quantiles ----
  
  qr <- quantile(subdatas$"rin", probs = seq(0, 1, 0.25))[c(2, 4)]
  qr <- qr
  names(qr) <- c("Q75R", "Q25R")
  
  qd <- quantile(subdatas$"din", probs = seq(0, 1, 0.25))[c(2, 4)]
  qd <- qd
  names(qd) <- c("Q25D", "Q75D")
  
  
  ## Create Other Category ----
  
  subdatas$"dr_class"[!(subdatas$"dr_class" %in% classes)] <- "Others"
  
  subdatas$"dr_class" <- factor(subdatas$"dr_class", levels = c(classes, "Others"))
  
  
  ## Add Scatterplot ----
  
  gplot <- ggplot(subdatas, aes(x = 1 - rin, y = din, col = dr_class))+
    
    geom_point(size = 3) + 
    
    scale_x_continuous(trans = "log10") + 
    
    scale_color_manual(values = c(color_classes, "#CCCCCC"))
  
  
  ## Get Plot Area Coordinates ----
  
  x_lims <- ggplot_build(gplot)$layout$panel_params[[1]]$x.range
  lag_x  <- diff(x_lims) * 0.04
  y_lims <- ggplot_build(gplot)$layout$panel_params[[1]]$y.range
  lag_y  <- diff(y_lims) * 0.04
  
  
  ## Add Rules ----
  
  gplot <- gplot +
    
    geom_vline(aes(xintercept = 1 - qr["Q75R"]), colour = "white", size = 2) + 
    geom_vline(aes(xintercept = 1 - qr["Q25R"]), colour = "white", size = 2) + 
    
    geom_hline(aes(yintercept = qd["Q75D"]), colour = "white", size = 2) + 
    geom_hline(aes(yintercept = qd["Q25D"]), colour = "white", size = 2) +
  
    geom_vline(aes(xintercept = 1 - qr["Q75R"]), colour = dark_grey) + 
    geom_vline(aes(xintercept = 1 - qr["Q25R"]), colour = dark_grey) + 
    
    geom_hline(aes(yintercept = qd["Q75D"]), colour = dark_grey) + 
    geom_hline(aes(yintercept = qd["Q25D"]), colour = dark_grey) 
  

  ## Change Theme ----
  
  gplot <- gplot +
    
    labs(x = "Restrictiveness", y = "Distinctiveness") +
    
    gg_theme +
    
    theme(legend.position = "None")

  
  ## Add DR Class Legend (left panel) ----
  
  if (taxa == "mammals") {
    
    segs <- data.frame(
      x     = rep(0.12, 3),
      xend  = rep(0.20, 3),
      y     = seq(0.1, 0.01, by = -0.035),
      label = classes_labs
    )
    
    gplot <- gplot +
      
      geom_rect(aes(xmin = segs[1, "x"] - 0.03, 
                    xmax = segs[1, "xend"] + 1.05, 
                    ymin = 0.00, 
                    ymax = 0.13),
                fill = "white", color = "white", size = 0.25)
    
    for (i in 1:nrow(segs)) {
      
      gplot <- gplot +
        
        annotate(geom = "segment", x = segs[i, "x"], y = segs[i, "y"],
                 xend = segs[i, "xend"], yend = segs[i, "y"], size = 1,
                 color = color_classes[i]) +
        
        annotate(geom = "text", x = segs[i, "xend"] + 0.05, y = segs[i, "y"], 
                 label = segs[i, "label"], hjust = "left", size = 8.5,
                 family  = par_family, colour = dark_grey, fontface = 1)
    }
  }
  
  
  ## Add Sub-Plot Label ----
  
  label <- ifelse(taxa == "mammals", "a", "b")
  
  gplot <- gplot +
    
    annotate(geom = "text", x = 0.00002, y = y_lims[2] - lag_y, 
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


options(warn = owarn)

## Message ----

usethis::ui_done(
  paste(
    usethis::ui_value(figname), 
    "successfully exported in figures/"
  )
)
