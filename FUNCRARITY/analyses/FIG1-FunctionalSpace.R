################################################################################
###                                                                          ###
###                        FIGURE 1 (PCoA SPACES)                            ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### AUTHORS : The Three Nico                                                 ###
### DATE    : 2019/08/05                                                     ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### > sessionInfo()                                                          ###
###                                                                          ###
### R version 3.5.3 (2019-03-11) -- "Great Truth"                            ###
### Platform: x86_64-apple-darwin18.2.0 (64-bit)                             ###
### Running under: macOS Mojave 10.14.4                                      ###
###                                                                          ###
### locale:                                                                  ###
### [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8        ###
###                                                                          ###
### attached base packages:                                                  ###
### [1] grid   stats  graphics  grDevices utils  datasets  methods   base    ###
###                                                                          ###
### other attached packages:                                                 ###
### [1] gridExtra_2.3  cowplot_1.0.0  png_0.1-7  ggtree_1.14.6               ###
### [5] ggplot2_3.1.0  treeio_1.6.2   ape_5.3                                ###
###                                                                          ###
################################################################################



rm(list = ls())



#'  -------------------------------------------------------------------------   @Parameters


root      <- "~/Desktop/FUNCRARITY/"
source(file.path(root, "graphsParameters.R"))

filename <- "Figure_1"

n        <- 1                          # ID of first plot
plots    <- list()                     # Subplots storage



#'  -------------------------------------------------------------------------   @LoadAddings


library(ggplot2)
library(png)
library(cowplot)
library(grid)
library(gridExtra)

addings <- list.files(path = path_R, pattern = "\\.R$", full.names = TRUE)
for (i in 1:length(addings)) { source(addings[i]) }



for (taxa in taxas) {



#' ---------------------------------------------------------------------------- @ImportData


  species_pcoa <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_pcoa.RData"))
    )
  )

  species_dr   <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_dr.RData"))
    )
  )

  species_fr     <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_funrar.RData"))
    )
  )

  species_silh  <- readPNG(
    source = file.path(path_data, taxa, paste0(taxa, "_silhouette.png"))
  )



#' ---------------------------------------------------------------------------- @SelectSpecies


  species_pcoa <- species_pcoa$vectors[
    rownames(species_pcoa$vectors) %in% rownames(species_dr), ]



#' ---------------------------------------------------------------------------- @SelectPCoAAxis


  data_all <- data.frame(
    x = jitter(species_pcoa[ , pcoa_axes[[taxa]][1]], jitter_val),
    y = jitter(species_pcoa[ , pcoa_axes[[taxa]][2]], jitter_val)
  )

  data_all <- merge(data_all, species_dr, by = "row.names")
  rownames(data_all) <- data_all[ , "Row.names"]

  data_all <- data_all[ , -1]

  data_all <- merge(data_all, species_fr$FR, by = "row.names", all.x = TRUE, all.y = FALSE)



#' ---------------------------------------------------------------------------- @PlotSpecies


  gplot <- ggplot(
    data     = data_all,
    mapping  = aes(x, y, fill = Din)
  ) +

  geom_point(
    shape   = 21,
    size    = 3,
    colour  = "transparent"
  ) +

  scale_fill_gradientn(colours = color_distinctiveness)



#' ---------------------------------------------------------------------------- @AddHulls


  for (classe in classes) {

    coords <- data_all[data_all[ , "DR_class"] == classe, ]
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
      # fill    = paste0(color_classes[classe], "22")
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
      label = c("Common", "Average", "Rare")
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
    image  = species_silh,
    x      = ifelse(taxa == "birds", -0.32, -0.30),
    y      = ifelse(taxa == "birds",  0.40,  0.40),
    scale  = ifelse(taxa == "mammals",  0.15,  0.12)
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
  filename  = file.path(path_figs, paste0(filename, ".png")),
  plot      = grobs,
  width     = 24,
  height    = 12,
  units     = "in",
  dpi       = 600
)
