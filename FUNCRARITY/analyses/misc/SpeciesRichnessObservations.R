library(raster)


root        <- "/Users/nicolascasajus/OneDrive/OneDrive - Fondation Biodiversité/MySpace/GROUPS/FREE/01-Loiseau/RALLL/FUNCRARITY/"
path_data   <- paste0(root, "data")
path_biomod <- "/Volumes/ROBERTO/free-sdm/"

taxas       <- c("mammals", "birds")
horizons    <- c("1979-2013", "2041-2060", "2061-2080")

study <- raster(file.path(path_data, "grid_area_with_cell_ids.tif"))



for (taxa in taxas) {

  species_taxo   <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_taxoinfos.RData"))
    )
  )

  projs <- list.files(
    path        = path_biomod,
    recursive   = TRUE,
    pattern     = paste0("^proj.+_1979-2013_bins$"),
    full.names  = TRUE
  )

  species <- strsplit(projs, "_")
  species <- unlist(
    lapply(
      species,
      function(x) {
        x[length(x)-2]
      }
    )
  )

  species <- species[which(species %in% species_taxo[ , "ID"])]


  spocc <- get(
    load(
      file.path(
        path_data,
        taxa,
        paste0("occurrences_", taxa)
      )
    )
  )

  spocc <- spocc[which(names(spocc) %in% species)]

  for (i in 1:length(spocc)) {

    cat("    Species:", i, "\r")

    cells <- data.frame(
      cell_id     = spocc[[i]],
      occurrence  = 1
    )


    ### Convert Sp Occurrence to SPDF -------------------

    pos <- which(!is.na(study[]))
    cell_id <- study[][pos]

    xy <- xyFromCell(study, pos)

    grille <- data.frame(
      cell_id = cell_id,
      x       = xy[ , 1],
      y       = xy[ , 2]
    )

    pts <- merge(cells, grille, by = "cell_id", all = FALSE)

    grid <- study
    pos <- which(!is.na(grid[]))
    grid[][pos] <- 0
    grid[][cellFromXY(grid, pts[ , c("x", "y")])] <- 1

    if (i == 1) {

      rich <- grid

    } else {

      rich[] <- rich[] + grid[]
    }
  }
  cat("\n")

  save(rich, file = paste0(path_data, "/", taxa, "/biomod-richness_", taxa, "_observations.RData"))
}
