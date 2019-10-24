rm(list = ls())



#'  -------------------------------------------------------------------------   @Parameters


root      <- "/Users/nicolascasajus/OneDrive/OneDrive - Fondation Biodiversité/MySpace/GROUPS/FREE/01-Loiseau/RALLL/FUNCRARITY/"
source(file.path(root, "graphsParameters.R"))

horizon  <- "2041-2060"



#'  -------------------------------------------------------------------------   @LoadAddings


library(vioplot)
library(png)

addings <- list.files(path = path_R, pattern = "\\.R$", full.names = TRUE)
for (i in 1:length(addings)) { source(addings[i]) }




#'  -------------------------------------------------------------------------   @ImportSpeciesSilhouettes


icons <- list()

for (taxa in taxas) {

  icons[[taxa]] <- readPNG(
      source = file.path(path_data, taxa, paste0(taxa, "_silhouette.png"))
  )
}



#'  -------------------------------------------------------------------------   @ImportData


threats <- get(
  load(
    file = file.path(path_data, "threats_nocc.RData")
  )
)

for (i in 2:ncol(threats)) { threats[ , i] <- as.character(threats[ , i]) }


climate <- get(
  load(
    file = file.path(path_data, "threats_cc.RData")
  )
)

for (i in 2:ncol(climate)) { climate[ , i] <- as.character(climate[ , i]) }

climate <- climate[climate[ , "Horizon"] == horizon, -ncol(climate)]


iucn <- get(
  load(
    file = file.path(path_data, "iucn_status.RData")
  )
)


#'  -------------------------------------------------------------------------   @PrepareData


datas <- rbind(threats, climate)

mat <- data.frame()

for (taxa in taxas) {

  tab <- datas[datas[ , "Taxa"] == taxa, ]

  for (threat in threats_vars) {

    dat <- tab[tab[ , "Threats"] == threat , ]

    mod     <- aov(Value ~ DR_class, data = dat)
    posthoc <- TukeyHSD(x = mod, which = "DR_class", conf.level = 0.95)
    posthoc <- posthoc$DR_class

    vars <- strsplit(rownames(posthoc), "-")
    var1 <- unlist(lapply(vars, function(x) x[1]))
    var2 <- unlist(lapply(vars, function(x) x[2]))

    tmp <- data.frame(
      taxa       = taxa,
      threat     = threat,
      var1       = var1,
      var2       = var2,
      pval       = posthoc[ , "p adj"],
      row.names  = NULL,
      stringsAsFactors = FALSE
    )

    mat <- rbind(mat, tmp)
  }
}

dat <- data.frame()
for (i in 1:nrow(mat)) {

  tmp <- rbind(mat[i, ], mat[i, ])
  dat <- rbind(dat, tmp)
}

for (i in seq(2, nrow(dat), 2)) {

  xx <- dat[i, "var1"]
  yy <- dat[i, "var2"]
  dat[i, "var1"] <- yy
  dat[i, "var2"] <- xx
}


mat <- data.frame()

for (taxa in taxas) {

  tab <- dat[dat[ , "taxa"] == taxa, ]

  for (threat in threats_vars) {

    kat <- tab[tab[ , "threat"] == threat , ]

    for (classe in classes) {

      tmp <- kat[kat[ , "var1"] == classe , ]

      xxx <- NULL

      for (i in 1:nrow(tmp)) {

        if (tmp[i, "pval"] < 0.05) {

          xxx <- c(xxx, tmp[i, "var2"])
        }
      }

      xxx <- gsub("D75R75", "a", xxx)
      xxx <- gsub("AVG", "b", xxx)
      xxx <- gsub("D25R25", "c", xxx)
      xxx <- paste0(sort(xxx), collapse = ",")

      tmp <- data.frame(
        taxa       = taxa,
        threat     = threat,
        DR_class   = classe,
        signif     = xxx,
        row.names  = NULL,
        stringsAsFactors = FALSE
      )

      mat <- rbind(mat, tmp)
    }
  }
}


save(mat, file = paste0(path_data, "/violin_pvalues_", horizon, ".RData"))
