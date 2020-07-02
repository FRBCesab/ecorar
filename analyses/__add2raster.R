ras   <- raster("data/grid_area_with_cell_ids.tif")
nulls <- get(load("data/SES_final.RData"))
nulls$cell <- as.numeric(as.character(nulls$cell))

pos <- which(ras[] %in% nulls$cell)

cells <- data.frame(pos = pos, cell = ras[][pos])
nulls <- merge(cells, nulls, by = "cell", all = TRUE)

birds_null = ras
birds_null[][nulls$pos] <- nulls$"SES_logmin.birds"

mammals_null = ras
mammals_null[][nulls$pos] <- nulls$"SES_logmin.mammals"
