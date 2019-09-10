
#'  -------------------------------------------------------------------------   @DefinePaths


root      <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/free-phylo/"
path_data <- paste0(root, "data/")



#'  -------------------------------------------------------------------------   @LoadAddings


library(ape)
library(funrar)
library(picante)
library(ggplot2)
library(gridExtra)




#'  -------------------------------------------------------------------------   @Parameters


n      <- 1                                         # ID of first plot
plots  <- list()                                    # Plot storage
vars   <- c("OriDin", "EDin")                       # Variables to plot
labels <- c("Originality", "Distinctiveness")       # Variable names

taxas  <- c("birds", "mammals")


for (taxa in taxas) {



#' ---------------------------------------------------------------------------- @ImportData


  species_ids    <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "ID.RData"))
    )
  )

  species_traits <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "trait.RData"))
    )
  )

  species_occs   <- get(
    load(
      file = file.path(path_data, taxa, paste0("occ_", taxa, "_list.RData"))
    )
  )

  species_fr     <- get(
    load(
      file = file.path(path_data, taxa, paste0("FR_", taxa, ".RData"))
    )
  )

  species_taxo   <- get(
    load(
      file = file.path(path_data, taxa, paste0("taxaInfo_", taxa, ".RData"))
    )
  )


  if (taxa == "mammals") {

    species_phylo <- get(
      load(
        file = file.path(path_data, taxa, paste0(taxa, "Phy.rdata"))
      )
    )

  }

  if (taxa == "birds") {

    species_phylo <- ape::read.tree(
      file = file.path(path_data, taxa, paste0(taxa, "_zillaHackett1_1.tre"))
    )
  }



#' ---------------------------------------------------------------------------- @SelectSpecies


  species_ids  <- species_ids[species_ids[ , "ID"] %in% rownames(species_traits), ]
  species_taxo <- species_taxo[species_taxo[ , "ID"] %in% species_ids[ , "ID"], ]

  species_set <- ape::drop.tip(
    phy = species_phylo,
    tip = species_phylo$tip.label[
      !is.element(
        el   = species_phylo$tip.label,
        set  = as.character(gsub(" ", "_", species_ids[ , "Name"]))
      )
    ]
  )



#' ---------------------------------------------------------------------------- @ComputePhylogeneticDistance


  species_distphyl <- cophenetic(species_set)
  species_distphyl <- sqrt(species_distphyl / max(species_distphyl))



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
    (species_fr$FR$Rin > q_r75) &
    (species_fr$FR$Rin < q_r25)
  ] <- "AVG"

  dr_classes <- data.frame(
    DR_class   = species_fr$FR[ , "DR_class"],
    row.names  = rownames(species_fr$FR)
  )



#' ---------------------------------------------------------------------------- @EvolutionnaryDistinctiveness


  sim_commu             <- matrix(1, nrow = 1, ncol = ncol(species_distphyl))
  colnames(sim_commu)   <- colnames(species_distphyl)

  species_edi <- t(funrar::distinctiveness(sim_commu, species_distphyl))
  species_edi <- data.frame(
    Escaped    = rownames(species_edi),
    EDi        = species_edi[ , 1],
    row.names  = NULL
  )


  species_ids <- as.data.frame(species_ids)
  species_ids[ , "Escaped"] <- as.character(gsub(" ", "_", species_ids[ , "Name"]))

  species_edi <- merge(
    x   = species_edi,
    y   = species_ids,
    by  = "Escaped"
  )

  species_edi$EDin      <- (species_edi[ , "EDi"] - min(species_edi[ , "EDi"])) /
                           max(species_edi[ , "EDi"] - min(species_edi[ , "EDi"]))

  rownames(species_edi) <- species_edi[ , "ID"]



#' ---------------------------------------------------------------------------- @EvolutionnaryOriginality


  ty <- picante::evol.distinct(
    tree = species_set,
    type = "equal.splits"
  )

  ty <- merge(x = ty, y = species_ids, by.x = "Species", by.y = "Escaped")

  ty[ , "tyin"] <- (ty[ , "w"] - min(ty[ , "w"])) / max(ty[ , "w"] - min(ty[ , "w"]))
  ty <- ty[ , c("ID", "w", "tyin")]



#' ---------------------------------------------------------------------------- @MergeData


  species_edi <- merge(
    x   = species_edi,
    y   = ty,
    by  = "ID"
  )
  species_edi <- species_edi[ , c("ID", "EDi", "EDin", "w", "tyin")]
  colnames(species_edi) <- c("ID", "EDi", "EDin", "OriDi", "OriDin")


  species_fr$FR[ , "ID"] <- rownames(species_fr$FR)
  species_fr$FR <- merge(
    x   = species_fr$FR,
    y   = species_edi,
    by  = "ID"
  )

  rownames(species_fr$FR) <- species_fr$FR[ , "ID"]
  species_fr$FR           <- species_fr$FR[ , -1]



#' ---------------------------------------------------------------------------- @BoxPlots


  dr_classes <- na.omit(species_fr$FR)


  for (j in 1:length(vars)) {

    plots[[n]] <- ggplot(
      data     = dr_classes,
      mapping  = aes_string(
        x     = "DR_class",
        y     = vars[j],
        fill  = "DR_class"
      )
    ) +

    geom_boxplot() +

    guides(fill = FALSE) +

    geom_jitter(
      width  = 0.1,
      size   = 0.5,
      color  = "darkgrey"
    ) +

    scale_y_continuous(limits = c(0, 1)) +

    geom_hline(
      yintercept  = mean(dr_classes[ , vars[j]], na.rm = TRUE),
      col         = "red",
      linetype    = "dashed"
    ) +

    labs(
      x = "DR class",
      y = paste("Evolutionary", labels[j])) +

    theme_bw()

    n <- n + 1
  }
}


do.call(gridExtra::grid.arrange, c(plots, ncol = 2))
