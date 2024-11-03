if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
}

library(SummarizedExperiment)
library(GEOquery)


# Llegir les dades
dades <- read.csv("human_cachexia.csv")

# Preparar dades per SummarizedExperiment
row_data <- dades[, c("Patient.ID", "Muscle.loss")]
assay_data <- as.matrix(dades[, -c(1, 2)])  # Excloem les dues primeres columnes, ja estan a row_data
col_data <- data.frame(Compound = colnames(assay_data)) # Metadades de les columnes

# Creació de l'objecte SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = assay_data),
  rowData = row_data,
  colData = col_data
)


# Exploració de les dades
summary(rowData(se)$`Muscle.loss`)  # Comptar quants pacients són cachectics o no
summary(assays(se)$counts)          # Resum de les concentracions de compostos

