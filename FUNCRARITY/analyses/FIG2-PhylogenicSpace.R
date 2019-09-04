################################################################################
###                                                                          ###
###                   FIGURE 2 (PHYLOGENIC TREES)                            ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### AUTHORS : The Three Nico                                                 ###
### DATE    : 2019/07/30                                                     ###
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


root      <- "/Users/nicolascasajus/OneDrive/OneDrive - Fondation Biodiversité/MySpace/GROUPS/FREE/01-Loiseau/RALLL/FUNCRARITY/"
source(file.path(root, "graphsParameters.R"))

filename <- "Figure_2"

n        <- 1                          # ID of first plot
n_lines  <- 9                          # Number of orders per column (in legend)
plots    <- list()                     # Subplots storage



#'  -------------------------------------------------------------------------   @LoadAddings


library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(png)
library(cowplot)
library(grid)
library(gridExtra)



for (taxa in taxas) {



#' ---------------------------------------------------------------------------- @ImportData


  species_ids    <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_ids.RData"))
    )
  )

  species_traits <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_traits.RData"))
    )
  )

  species_fr     <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_funrar.RData"))
    )
  )

  species_taxo   <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_taxoinfos.RData"))
    )
  )

  species_phylo <- ape::read.tree(
    file = file.path(path_data, taxa, paste0(taxa, "_phylo.tre"))
  )

  species_silh  <- readPNG(
    source = file.path(path_data, taxa, paste0(taxa, "_silhouette.png"))
  )

  d_d75r75      <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_d_d75r75.RData"))
    )
  )

  d_avg         <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_d_avg.RData"))
    )
  )

  d_d25r25      <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_d_d25r25.RData"))
    )
  )

  species_d     <- data.frame(
    rbind(d_d75r75, d_avg, d_d25r25),
    DR_class = c(
      rep("D75R75", nrow(d_d75r75)),
      rep("AVG",    nrow(d_avg)),
      rep("D25R25", nrow(d_d25r25))
    )
  )




#' ---------------------------------------------------------------------------- @SelectSpecies


  species_ids  <- species_ids[species_ids[ , "ID"] %in% rownames(species_traits), ]
  species_taxo <- species_taxo[species_taxo[ , "ID"] %in% species_ids[ , "ID"], ]

  species_ids  <- species_ids[species_ids[ , "ID"] %in% rownames(species_fr$FR[!is.na(species_fr$FR[ , "Din"]), ]), ]
  species_taxo <- species_taxo[species_taxo[ , "ID"] %in% rownames(species_fr$FR[!is.na(species_fr$FR[ , "Din"]), ]), ]

  species_set  <- ape::drop.tip(
    phy = species_phylo,
    tip = species_phylo$tip.label[
      !is.element(
        el   = species_phylo$tip.label,
        set  = as.character(gsub(" ", "_", species_ids[ , "Name"]))
      )
    ]
  )

  species_ids   <- species_ids[species_ids[ , "Name"] %in% gsub("_", " ", species_set$tip.label), ]
  species_taxo  <- species_taxo[species_taxo[ , "Name"] %in% gsub("_", " ", species_set$tip.label), ]
  species_fr$FR <- species_fr$FR[rownames(species_fr$FR) %in% species_ids[ , "ID"], ]



#' ---------------------------------------------------------------------------- @AddDRClasses


  q_d75 <- species_fr$Q$Q75_D
  q_d25 <- species_fr$Q$Q25_D
  q_r75 <- species_fr$Q$Q75_R
  q_r25 <- species_fr$Q$Q25_R

  species_fr$FR$DR_class[
    (species_fr$FR$Din < q_d25) &
    (species_fr$FR$Rin < q_r75)
  ] <- "D25R25"

  species_fr$FR$DR_class[
    (species_fr$FR$Din > q_d75) &
    (species_fr$FR$Rin > q_r25)
  ] <- "D75R75"

  species_fr$FR$DR_class[
    (species_fr$FR$Din < q_d25) &
    (species_fr$FR$Rin > q_r25)
  ] <- "D25R75"

  species_fr$FR$DR_class[
    (species_fr$FR$Din > q_d75) &
    (species_fr$FR$Rin < q_r75)
  ] <- "D75R25"

  species_fr$FR$DR_class[
    (species_fr$FR$Din > q_d25) &
    (species_fr$FR$Din < q_d75) &
    (species_fr$FR$Rin > q_r25) &
    (species_fr$FR$Rin < q_r75)
  ] <- "AVG"



#' ---------------------------------------------------------------------------- @AddData2PhyloObj


  species_taxo[ , "Name"] <- as.character(species_taxo[ , "Name"])
  rownames(species_taxo)  <- species_taxo[ , "ID"]

  groupInfo  <- tapply(
    X      = species_taxo[ , "Name"],
    INDEX  = species_taxo[ , "order"],
    function(x) { return(x) }
  )

  ordre  <- sapply(
    gsub("_", " ", species_set$tip.label),
    function(x, y) { return(which(y == x)) },
    y = species_taxo[ , "Name"]
  )

  ordres <- as.character(species_taxo[ordre, "order"])
  ordres <- unique(ordres)

  groupInfo <- groupInfo[ordres]

  groupInfo <- lapply(groupInfo, function(x) gsub(" ", "_", x))

  new_phylo <- groupOTU(species_set, groupInfo)


  species_fr$FR <- merge(
    x    = species_fr$FR,
    y    = species_taxo,
    by   = "row.names",
    all  = TRUE
  )

  species_fr$FR$Name <- gsub(" ", "_", species_fr$FR$Name)

  dat <- species_fr$FR[ , c("Name", "Din", "DR_class")]
  dat <- dat[match(gsub(" ", "_", names(ordre)), as.character(dat[ , "Name"])), ]
  dat$node <- 1:nrow(dat)

  dat["D25R25"] <- ifelse(dat[ , "DR_class"] == "D25R25", 1, NA)
  dat["D75R75"] <- ifelse(dat[ , "DR_class"] == "D75R75", 1, NA)
  dat["AVG"]    <- ifelse(dat[ , "DR_class"] == "AVG", 1, NA)

  dat <- tidytree::as_tibble(dat)

  new_phylo <- tidytree::as_tibble(new_phylo)
  new_phylo <- treeio::full_join(new_phylo, dat, by = "node")
  new_phylo <- tidytree::as.treedata(new_phylo)



#' ---------------------------------------------------------------------------- @CreateGGTree


  tree_plot <- ggtree(
    tr          = new_phylo,
    mapping     = aes(color = Din),
    layout      = "circular",
    ladderize   = FALSE,
    right       = TRUE,
    size        = 0.1
  ) +

  theme(
    plot.margin = unit(rep(ifelse(taxa == "birds", -3.00, -3.90), 4), "cm")
  ) +

  scale_colour_gradientn(colours = color_distinctiveness)



#' ---------------------------------------------------------------------------- @GetSegmentsCoords


  tree_dt <- as.data.frame(tree_plot$data)
  tree_dt <- tree_dt[tree_dt[ , "isTip"] == TRUE, ]
  tree_dt <- tree_dt[order(tree_dt[ , "y"]), ]

  groups <- unique(as.character(tree_dt[ , "group"]))
  n_groups <- length(groups)

  coord_groups <- data.frame()

  if (taxa == "birds") {

    space <- c(20, 15, 10)

  } else {

    space <- c(42, 32, 22)
  }

  space <- rep(space, round(n_groups / length(space)))

  for (i in 1:n_groups) {

    dat <- as.data.frame(matrix(nrow = 1, ncol = 10))
    colnames(dat) <- c(
      "group", "id_gr", "y1", "y2", "angle", "angle_adj", "n", "y_mid", "h_just", "x"
    )

    tmp <- tree_dt[tree_dt[ , "group"] == groups[i], ]

    dat["group"] <- groups[i]
    dat["id_gr"] <- i
    dat["y1"]    <- min(tmp[ , "y"])
    dat["y2"]    <- max(tmp[ , "y"])
    dat["angle"] <- mean(tmp[ , "angle"])
    dat["n"]     <- nrow(tmp)
    dat["y_mid"] <- mean(c(max(tmp[ , "y"]), min(tmp[ , "y"])))

    if (dat["n"] == 1) {

      dat["y1"] <- dat["y1"] - 0.1
      dat["y2"] <- dat["y2"] + 0.1
    }

    dat["angle_adj"] <- dat["angle"]

    if (dat["angle"] >= 90 && dat["angle"] <= 180) {

      dat["angle_adj"] <- dat["angle"] + 180

    } else {

      if (dat["angle"] > 180 && dat["angle"] <= 270) {

        dat["angle_adj"] <- dat["angle"] - 180
      }
    }

    dat["h_just"] <- ifelse(dat["angle"] >= 90 && dat["angle"] <= 270, 1L, 0L)

    dat["x"] <- max(tree_dt["x"]) + space[i]

    coord_groups <- rbind(coord_groups, dat)
  }



#' ---------------------------------------------------------------------------- @AddSegments


  tree_plot <- tree_plot +

  geom_segment(
    data     = coord_groups,
    mapping  = aes(
      x     = x,
      y     = y1,
      xend  = x,
      yend  = y2
    ),
    color    = light_grey,
    lineend  = "butt",
    size     = 0.5
  )



#' ---------------------------------------------------------------------------- @AddSegmentsLabels


  tree_plot <- tree_plot +

  geom_text(
    data     = coord_groups,
    mapping  = aes(
      x      = x,
      y      = y_mid,
      hjust  = h_just,
      label  = id_gr
    ),
    vjust    = 0.45,
    size     = 4.0,
    nudge_x  = 1.25,
    color    = dark_grey
  )



#' ---------------------------------------------------------------------------- @AddD75R75Points


  cols <- ifelse(is.na(tree_dt[ , "D75R75"]), NA, paste0(color_rare, alpha))

  tree_dt[ , "x"] <- max(tree_dt[ , "x"]) + ifelse(taxa == "birds", 1.5, 2.5)

  tree_plot <- tree_plot +

  geom_point(
    data     = tree_dt,
    mapping  = aes(
      x      = x,
      y      = y,
      color  = D75R75
    ),
    fill     = cols,
    color    = "transparent",
    shape    = 21,
    size     = 3
  )



#' ---------------------------------------------------------------------------- @AddAVGPoints


  cols <- ifelse(is.na(tree_dt[ , "AVG"]), NA, paste0(color_avg, alpha))

  tree_dt[ , "x"] <- max(tree_dt[ , "x"]) + ifelse(taxa == "birds", 2.5, 4.0)

  tree_plot <- tree_plot +

  geom_point(
    data     = tree_dt,
    mapping  = aes(
      x      = x,
      y      = y,
      color  = AVG
    ),
    fill     = cols,
    color    = "transparent",
    shape    = 21,
    size     = 3
  )



#' ---------------------------------------------------------------------------- @AddD25R25Points


  cols <- ifelse(is.na(tree_dt[ , "D25R25"]), NA, paste0(color_common, alpha))

  tree_dt[ , "x"] <- max(tree_dt[ , "x"]) + ifelse(taxa == "birds", 2.5, 4.0)

  tree_plot <- tree_plot +

  geom_point(
    data     = tree_dt,
    mapping  = aes(
      x      = x,
      y      = y,
      color  = D25R25
    ),
    fill     = cols,
    color    = "transparent",
    shape    = 21,
    size     = 3
  )



#' ---------------------------------------------------------------------------- @AddSilhouette


  tree_plot <- cowplot::ggdraw(tree_plot) +

  cowplot::draw_image(
    image  = species_silh,
    x      = ifelse(taxa == "mammals", -0.37, -0.38),
    y      = ifelse(taxa == "mammals",  0.42,  0.42),
    scale  = ifelse(taxa == "mammals",  0.15,  0.12)
  )



#' ---------------------------------------------------------------------------- @AddHistogram


  hist_plot <- ggplot(
    data     = species_d,
    mapping  = aes(
      x      = estimated_D,
      color  = DR_class,
      fill   = DR_class
    )
  ) +

  geom_density(adjust = 1.5) +

  scale_x_continuous(limits = c(0, 1)) +

  scale_color_manual(
    values = c(color_avg, color_common, color_rare)
  ) +

  scale_fill_manual(
    values = paste0(c(color_avg, color_common, color_rare), alpha)
  ) +

  theme_light() +

  theme(
    axis.title       = element_blank(),
    legend.position  = "None"
  ) +

 annotate(
   geom    = "text",
   x       = 0.25,
   y       = 47.5,
   label   = "bold(\"Density of D\")",
   color   = dark_grey,
   size    = 5,
   family  = "serif",
   parse = TRUE
 )

  hist_plot <- ggplotGrob(hist_plot)

  tree_plot <- tree_plot +

  annotation_custom(
    grob = hist_plot,
    xmin = 0.40,
    xmax = 0.60,
    ymin = 0.40,
    ymax = 0.60
  )


  plots[[n]] <- tree_plot

  n <- n + 1



#' ---------------------------------------------------------------------------- @AddOrderLegend


  tree_legend <- ggplot() +

  theme_bw() +

  theme(
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    line              = element_blank(),
    text              = element_blank(),
    title             = element_blank()
  )

  n_columns <- ceiling(nrow(coord_groups) / n_lines)

  xx    <- ifelse(taxa == "mammals", 11, 1)
  pos   <- 0

  for (j in 1:n_columns) {

    for (yy in n_lines:1) {

      pos <- pos + 1

      if (pos <= nrow(coord_groups)) {

        texte <- coord_groups[pos, "group"]
        substr(texte, 1, 1)            <- toupper(substr(texte, 1, 1))
        substr(texte, 2, nchar(texte)) <- tolower(substr(texte, 2, nchar(texte)))
        texte <- paste0(coord_groups[pos, "id_gr"], "  ", texte)

        tree_legend <- tree_legend +

        annotate(
          geom    = "text",
          x       = xx,
          y       = yy,
          size    = 5,
          label   = texte,
          hjust   = "left",
          color   = dark_grey,
          family  = "serif"
        )
      }
    }

    xx <- xx + ifelse(taxa == "mammals", 6.75, 5.80)
  }


  coords <- data.frame(
    xmin = ifelse(taxa == "mammals", 10, 0),
    xmax = 30,
    ymin = 0,
    ymax = n_lines + 1
  )

  tree_legend <- tree_legend +

  geom_rect(
    data     = coords,
    mapping  = aes(
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
    ),
    fill     = "transparent",
    color    = dark_grey,
    size     = 0.25
  ) +

  scale_x_continuous(limits = c(0, 30)) +

  scale_y_continuous(limits = c(0, n_lines + 1))



#' ---------------------------------------------------------------------------- @AddColorsLegend


  if (taxa == "mammals") {

    coords <- data.frame(
      x       = rep(1.0, 3),
      x_text  = rep(1.5, 3),
      y       = seq(n_lines - 1, n_lines - 3, by = -1),
      text    = c("Rare species", "Average", "Common species")
    )

    yctr <- 2.85
    xs <- seq(0.9, 6.1, length.out = length(color_distinctiveness) + 1)
    gradient <- data.frame(
      x1 = xs[-length(xs)],
      x2 = xs[-1],
      y1 = rep(yctr - .5, length(xs) - 1),
      y2 = rep(yctr + .5, length(xs) - 1)
    )

    tree_legend <- tree_legend +

    geom_point(
      data     = coords,
      mapping  = aes(
        x      = x,
        y      = y
      ),
      fill     = c(color_rare, color_avg, color_common),
      color    = "transparent",
      shape    = 21,
      size     = 4
    ) +

    geom_text(
      data     = coords,
      mapping  = aes(
        x      = x_text,
        y      = y,
        hjust  = -2,
        label  = text
      ),
      hjust    = "left",
      size     = 5.0,
      color    = dark_grey,
      family   = "serif"
    ) +

    geom_rect(
      mapping  = aes(
        xmin = 0, xmax = 7, ymin = 0, ymax = n_lines + 1
      ),
      fill     = "transparent",
      color    = dark_grey,
      size     = 0.25
    ) +

    geom_rect(
      data     = gradient,
      mapping  = aes(
        xmin = x1,
        xmax = x2,
        ymin = y1,
        ymax = y2
      ),
      fill     = color_distinctiveness
    ) +

    geom_rect(
      mapping  = aes(
        xmin = min(xs), xmax = max(xs), ymin = yctr - .5, ymax = yctr + .5
      ),
      fill     = "transparent",
      color    = dark_grey,
      size     = 0.25
    ) +

    annotate(
      geom    = "text",
      x       = xs[1],
      y       = yctr - 1,
      label   = "0",
      color   = dark_grey,
      size    = 4,
      family  = "serif"
    ) +

    annotate(
      geom    = "text",
      x       = xs[length(xs)],
      y       = yctr - 1,
      label   = "1",
      color   = dark_grey,
      size    = 4,
      family  = "serif"
    ) +

    annotate(
      geom    = "text",
      x       = 3.5,
      y       = yctr + 1.15,
      label   = "bold(Distinctiveness)",
      color   = dark_grey,
      size    = 5,
      family  = "serif",
      parse = TRUE
    )
  }


  plots[[n]] <- tree_legend

  n <- n + 1
}



#' ---------------------------------------------------------------------------- @ArrangeSubplots


mat <- matrix(
  data   = c(rep(c(rep(1, 4), rep(3, 4)), 4), rep(2, 4), rep(4, 4)),
  ncol   = 8,
  nrow   = 5,
  byrow  = TRUE
)

grobs <- gridExtra::arrangeGrob(
  plots[[1]], plots[[2]], plots[[3]], plots[[4]],
  layout_matrix = mat
)



#' ---------------------------------------------------------------------------- @ExportFigure


ggsave(
  filename  = file.path(path_figs, paste0(filename, ".png")),
  plot      = grobs,
  width     = 27,
  height    = 15,
  units     = "in",
  dpi       = 600
)



#' ---------------------------------------------------------------------------- @NOTES


# https://yulab-smu.github.io/treedata-book/index.html
# https://rawgit.com/valentinitnelav/valentinitnelav.github.io/master/assets/2018-01-07-ggtree/2018-01-07-ggtree.html
