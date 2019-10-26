#' --------------------------------------------------------------------------   @Header
#'
#' @title Figure 2 - Panel of Phylogeny Trees
#'
#' @description
#' This R script...
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#'         Nicolas LOISEAU, \email{nicolas.loiseau1@@gmail.com}
#'
#' @date 2019/08/05
#'
#' --------------------------------------------------------------------------   @VVVVVV




#'  -------------------------------------------------------------------------   @Parameters


n        <- 1                          # ID of first plot
n_lines  <- 9                          # Number of orders per column (in legend)
plots    <- list()                     # Subplots storage



for (taxa in taxas) {



  species_d <- estimated_d[estimated_d[ , "class"] == taxa, ]
  

#' ---------------------------------------------------------------------------- @ImportData


  species_phylo <- ape::read.tree(
    file = file.path(path_data, paste0(taxa, "_phylo.tre"))
  )



#' ---------------------------------------------------------------------------- @SelectSpecies


  species_set  <- ape::drop.tip(
    phy = species_phylo,
    tip = species_phylo$tip.label[
      !is.element(
        el   = species_phylo$tip.label,
        set  = gsub("\\s", "_", datas[ , "scientific_name"])
      )
    ]
  )

  species_info <- datas
  species_info[ , "grep_name"] <- gsub("\\s", "_", species_info[ , "scientific_name"])

  species_info <- species_info[species_info[ , "grep_name"] %in% species_set$tip.label, ]



#' ---------------------------------------------------------------------------- @SortSpeciesByOrder


  group_info  <- tapply(
    X     = species_info[ , "grep_name"],
    INDEX = species_info[ , "order"],
    function(x) { return(x) }
  )

  ordre <- sapply(
    species_set$tip.label,
    function(x, y) { return(which(y == x)) },
    y = species_info[ , "grep_name"]
  )

  ordres <- unique(species_info[ordre, "order"])

  group_info <- group_info[ordres]

  group_info <- lapply(group_info, function(x) x)

  new_phylo <- groupOTU(species_set, group_info)



#' ---------------------------------------------------------------------------- @AddData2PhyloObj


  dat <- species_info[ , c("grep_name", "din", "dr_class")]

  dat <- dat[match(names(ordre), dat[ , "grep_name"]), ]
  dat$node <- 1:nrow(dat)

  dat["D25R25"] <- ifelse(dat[ , "dr_class"] == "D25R25", 1, NA)
  dat["D75R75"] <- ifelse(dat[ , "dr_class"] == "D75R75", 1, NA)
  dat["AVG"]    <- ifelse(dat[ , "dr_class"] == "AVG", 1, NA)

  dat <- tidytree::as_tibble(dat)

  new_phylo <- tidytree::as_tibble(new_phylo)
  new_phylo <- treeio::full_join(new_phylo, dat, by = "node")
  new_phylo <- tidytree::as.treedata(new_phylo)



#' ---------------------------------------------------------------------------- @CreateGGTree


  tree_plot <- ggtree(
    tr          = new_phylo,
    mapping     = aes(color = din),
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
    image  = icons[[taxa]],
    x      = ifelse(taxa == "mammals", -0.37, -0.38),
    y      = ifelse(taxa == "mammals",  0.42,  0.42),
    scale  = ifelse(taxa == "mammals",  0.15,  0.12)
  )



#' ---------------------------------------------------------------------------- @AddHistogram


  hist_plot <- ggplot(
    data     = species_d,
    mapping  = aes(
      x      = estimated_d,
      color  = dr_class,
      fill   = dr_class
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
  filename  = file.path(path_figs, paste0(figname2, ".png")),
  plot      = grobs,
  width     = 27,
  height    = 15,
  units     = "in",
  dpi       = 600
)



#' ---------------------------------------------------------------------------- @NOTES


# https://yulab-smu.github.io/treedata-book/index.html
# https://rawgit.com/valentinitnelav/valentinitnelav.github.io/master/assets/2018-01-07-ggtree/2018-01-07-ggtree.html
