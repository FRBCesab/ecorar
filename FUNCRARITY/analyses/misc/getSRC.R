indir <- "/Volumes/ROBERTO/outputs#1/metrics/"


horizon <- "2061-2080"

fls <- list.files(path = indir, pattern = paste0("^src.+", horizon,"_bins"), full.names = FALSE)

# tab <- data.frame()

for (i in 1:length(fls)) {

  cat(i, "              \r")

  tmp <- get(load(paste0(indir, fls[i])))

  dat <- data.frame(
    species       = unlist(strsplit(fls[i], "_"))[2],
    horizon       = horizon,
    climateChange = tmp$"Compt.By.Models"[1, "SpeciesRangeChange"],
    stringsAsFactors = FALSE
  )

  tab <- rbind(tab, dat)

}


root      <- "/Users/nicolascasajus/OneDrive/OneDrive - Fondation Biodiversité/MySpace/GROUPS/FREE/01-Loiseau/FUNCRARITY/"
path_data <- paste0(root, "data")

species_taxo   <- rbind(
  data.frame(
    get(
      load(
        file = file.path(path_data, "mammals", "mammals_taxoinfos.RData")
      )
    ),
    taxa = "mammals"
  ),
  data.frame(
    get(
      load(
        file = file.path(path_data, "birds", "birds_taxoinfos.RData")
      )
    ),
    taxa = "birds"
  )
)

species_taxo <- species_taxo[ , c("taxa", "ID")]

for (i in 1:2) species_taxo[ , i] <- as.character(species_taxo[ , i])


src <- merge(species_taxo, tab, by.x = "ID", by.y = "species", all = TRUE)

pos <- which(!is.na(src[ , "climateChange"]))
if (length(pos) > 0) src <- src[pos, ]

pos <- which(src[ , "climateChange"] <= 500)
if (length(pos) > 0) src <- src[pos, ]


species_dr   <- rbind(
  data.frame(
    get(
      load(
        file = file.path(path_data, "mammals", "mammals_dr.RData")
      )
    ),
    taxa = "mammals"
  ),
  data.frame(
    get(
      load(
        file = file.path(path_data, "birds", "birds_dr.RData")
      )
    ),
    taxa = "birds"
  )
)

pos <- which(!is.na(species_dr[ , "DR_class"]))
if (length(pos) > 0) species_dr <- species_dr[pos, ]

pos <- which(species_dr[ , "DR_class"] %in% c("D75R75", "AVG", "D25R25"))
if (length(pos) > 0) species_dr <- species_dr[pos, ]

species_dr$ID <- rownames(species_dr)

src <- merge(src, species_dr, by = "ID", all = FALSE)


src2050 <- src[src[ , "horizon"] == "2041-2060", c("ID", "taxa.x", "DR_class", "climateChange")]
colnames(src2050) <- c("ID", "taxa", "class", "src2050")

src2080 <- src[src[ , "horizon"] == "2061-2080", c("ID", "climateChange")]
colnames(src2080) <- c("ID", "src2080")

src <- merge(src2050, src2080, by = "ID", all = FALSE)


occ <- read.delim(paste0(path_data, "/species_list_prevalence.txt"), header = TRUE)
occ <- occ[ , -1]
colnames(occ) <- c("ID", "occurrence")

src <- merge(src, occ, by = "ID", all = FALSE)

save(src, file = paste0(path_data, "/src_vs_occurrence.RData"))



# src <- src[ , c("climateChange", "DR_class", "taxa.x", "horizon")]
# src$Threats <- "Climate change"
#
# src <- src[ , c(1, 2, 5, 3, 4)]
# colnames(src) <- c("Value", "DR_class", "Threats", "Taxa", "Horizon")
#
# threats_cc <- src
#
# save(threats_cc, file = paste0(path_data, "threats_cc.RData"))
