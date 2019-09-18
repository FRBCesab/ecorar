library(raster)


root        <- "/Users/nicolascasajus/OneDrive/OneDrive - Fondation Biodiversité/MySpace/GROUPS/FREE/01-Loiseau/RALLL/FUNCRARITY/"
path_data   <- paste0(root, "data")
path_biomod <- "/Volumes/ROBERTO/free-sdm/"

taxas       <- c("mammals", "birds")
horizons    <- c("1979-2013", "2041-2060", "2061-2080")


for (taxa in taxas) {

  species_taxo   <- get(
    load(
      file = file.path(path_data, taxa, paste0(taxa, "_taxoinfos.RData"))
    )
  )

  for (horizon in horizons) {

    cat(">>>", taxa, "-", horizon, "\n")

    projs <- list.files(
      path        = path_biomod,
      recursive   = TRUE,
      pattern     = paste0("^proj.+_", horizon, "_bins$"),
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

    pos <- which(species %in% species_taxo[ , "ID"])
    projs <- projs[pos]



    for (i in 1:length(projs)) {

      cat("    Species:", i, "\r")

      tmp <- get(load(projs[i]))

      if (i == 1) {

        rich <- tmp

      } else {

        rich[] <- rich[] + tmp[]
      }
    }
    cat("\n")

    save(rich, file = paste0(path_data, "/", taxa, "/biomod-richness_", taxa, "_", horizon, ".RData"))
  }
}
