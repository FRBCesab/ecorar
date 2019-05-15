path_data <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/data/"           ### !!!
path_R    <- "~/OneDrive/OneDrive - Fondation Biodiversité/MySpace/RALLL/scripts/"  ### !!!

foot  <- raster(paste0(path_data, "HFP_2009.tif"))

load(paste0(path_data, "birds/occ_birds_list.RData"))
load(paste0(path_data, "birds/data_DR_birds.RData"))

data_DR_birds$DR_class <- as.character(data_DR_birds$DR_class)

pos <- which(data_DR_birds$DR_class == "D75R75")
d75r75 <- rownames(data_DR_birds)[pos]

pos <- which(data_DR_birds$DR_class == "D25R25")
d25r25 <- rownames(data_DR_birds)[pos]

pos <- which(data_DR_birds$DR_class == "AVG")
avg <- rownames(data_DR_birds)[pos]

cat("BIRDS - D75R75\n")

cells <- NULL
for (i in 1:length(d75r75)) {
  species <- d75r75[i]
  x <- lapply(occ_birds_list, function(x, species) ifelse(species %in% x, 1, 0), species)
  cells[i] <- mean(foot[][as.numeric(names(occ_birds_list)[which(x == 1)])], na.rm = TRUE)
}

hfps <- data.frame(
  Value     = cells,
  Threats   = "Mean HFP",
  DR_class  = "D75R75",
  Taxa      = "birds",
  stringsAsFactors = FALSE
)

cat("BIRDS - D25R25\n")

cells <- NULL
for (i in 1:length(d25r25)) {
  species <- d25r25[i]
  x <- lapply(occ_birds_list, function(x, species) ifelse(species %in% x, 1, 0), species)
  cells[i] <- mean(foot[][as.numeric(names(occ_birds_list)[which(x == 1)])], na.rm = TRUE)
}

subhfps <- data.frame(
  Value     = cells,
  Threats   = "Mean HFP",
  DR_class  = "D25R25",
  Taxa      = "birds",
  stringsAsFactors = FALSE
)
hfps <- rbind(hfps, subhfps)

cat("BIRDS - AVG\n")

cells <- NULL
for (i in 1:length(avg)) {
  species <- avg[i]
  x <- lapply(occ_birds_list, function(x, species) ifelse(species %in% x, 1, 0), species)
  cells[i] <- mean(foot[][as.numeric(names(occ_birds_list)[which(x == 1)])], na.rm = TRUE)
}

subhfps <- data.frame(
  Value     = cells,
  Threats   = "Mean HFP",
  DR_class  = "AVG",
  Taxa      = "birds",
  stringsAsFactors = FALSE
)
hfps <- rbind(hfps, subhfps)


####

load(paste0(path_data, "mammals/occ_mammals_list.RData"))
load(paste0(path_data, "mammals/data_DR_mammals.RData"))

data_DR_mammals$DR_class <- as.character(data_DR_mammals$DR_class)

pos <- which(data_DR_mammals$DR_class == "D75R75")
d75r75 <- rownames(data_DR_mammals)[pos]

pos <- which(data_DR_mammals$DR_class == "D25R25")
d25r25 <- rownames(data_DR_mammals)[pos]

pos <- which(data_DR_mammals$DR_class == "AVG")
avg <- rownames(data_DR_mammals)[pos]

cat("MAMMALS - D75R75\n")

cells <- NULL
for (i in 1:length(d75r75)) {
  species <- d75r75[i]
  x <- lapply(occ_mammals_list, function(x, species) ifelse(species %in% x, 1, 0), species)
  cells[i] <- mean(foot[][as.numeric(names(occ_mammals_list)[which(x == 1)])], na.rm = TRUE)
}

subhfps <- data.frame(
  Value     = cells,
  Threats   = "Mean HFP",
  DR_class  = "D75R75",
  Taxa      = "mammals",
  stringsAsFactors = FALSE
)
hfps <- rbind(hfps, subhfps)

cells <- NULL
for (i in 1:length(d25r25)) {
  species <- d25r25[i]
  x <- lapply(occ_mammals_list, function(x, species) ifelse(species %in% x, 1, 0), species)
  cells[i] <- mean(foot[][as.numeric(names(occ_mammals_list)[which(x == 1)])], na.rm = TRUE)
}

cat("MAMMALS - D25R25\n")

subhfps <- data.frame(
  Value     = cells,
  Threats   = "Mean HFP",
  DR_class  = "D25R25",
  Taxa      = "mammals",
  stringsAsFactors = FALSE
)
hfps <- rbind(hfps, subhfps)

cat("MAMMALS - AVG\n")

cells <- NULL
for (i in 1:length(avg)) {
  species <- avg[i]
  x <- lapply(occ_mammals_list, function(x, species) ifelse(species %in% x, 1, 0), species)
  cells[i] <- mean(foot[][as.numeric(names(occ_mammals_list)[which(x == 1)])], na.rm = TRUE)
}

subhfps <- data.frame(
  Value     = cells,
  Threats   = "Mean HFP",
  DR_class  = "AVG",
  Taxa      = "mammals",
  stringsAsFactors = FALSE
)
hfps <- rbind(hfps, subhfps)


save(hfps, file = "threat_hfp.RData")
