

#'  -------------------------------------------------------------------------   @DefinePaths


path_R    <- "R"
path_data <- "data"
path_figs <- "figures"



#'  -------------------------------------------------------------------------   @DefineVariables


taxas          <- c("mammals", "birds")

classes        <- c("D25R25", "AVG", "D75R75")
classes_labs   <- c("Common", "Average", "Rare")
vars_richness  <- c("TD_sp", "D75R75", "D25R25")

threats_vars    <- c("Protection Target", "meanHDI", "Human FootPrint", "Climate change")
threats_labels <- c("Target achievement (%)", "Mean HDI", "Human Footprint", "Climate change (%)")
iucn_status    <- c("NE", "LC", "TH")




#'  -------------------------------------------------------------------------   @DefineColors


color_rare   <- "#FF4500"             # ~ Red
color_avg    <- "#00AFBB"             # ~ Turquoise
color_common <- "#E7B800"             # ~ Orange

color_classes <- c(color_common, color_avg, color_rare)
names(color_classes) <- classes

alpha        <- "88"

light_grey   <- "#888888"
dark_grey    <- "#333333"

color_distinctiveness <- RColorBrewer::brewer.pal(name = "YlGnBu", n = 9)
color_distinctiveness <- colorRampPalette(color_distinctiveness)(255)

color_richness <- RColorBrewer::brewer.pal(name = "YlOrRd", n = 9)
color_richness <- colorRampPalette(color_richness)(255)
color_richness <- c("#AAAAAA", color_richness)

color_ocean <- "#95D8EB"

color_iucn_bg <- c("#AAAAAA", "#026666", "#C53131")
color_iucn_fg <- c("#333333", "#F7F7F7", "#F6C8C8")
names(color_iucn_bg) <- names(color_iucn_fg) <- iucn_status



#'  -------------------------------------------------------------------------   @OthersParameters


proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

pcoa_axes <- list()
pcoa_axes[["mammals"]] <- c(2, 3)
pcoa_axes[["birds"]] <- c(2, 4)

jitter_val <- 500
