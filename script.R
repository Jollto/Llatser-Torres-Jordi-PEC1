if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
}

library(SummarizedExperiment)
library(FactoMineR)
library(factoextra)

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
table(rowData(se)$Muscle.loss)  # Comptar quants pacients són cachectics o no
summary(assays(se)$counts)          # Resum de les concentracions de compostos


## Transformació logarítmica per evitar outliers
assays(se)$log_counts <- log(assays(se)$counts + 1)
summary(assays(se)$log_counts)

assays(se)$exploration_counts <- scale(assays(se)$log_counts)
summary(assays(se)$exploration_counts)

dades_pca <- princomp(assays(se)$exploration_counts)
fviz_pca_ind(dades_pca, 
             geom.ind = "point", 
             col.ind = rowData(se)$Muscle.loss, # Color per condició de pèrdua muscular
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Muscle loss") +
  labs(title = "PCA de les dades de concentració de compostos")



# Extreure les contribucions de les variables per a cada component principal
var_contrib <- get_pca_var(dades_pca)$contrib

# Seleccionar els 5 compostos amb més contribució per cada component
top_contrib_pc1 <- order(var_contrib[, "Dim.1"], decreasing = TRUE)[1:5]
top_contrib_pc2 <- order(var_contrib[, "Dim.2"], decreasing = TRUE)[1:5]
# Obtenir els noms dels compostos amb major contribució
top_vars <- unique(c(top_contrib_pc1, top_contrib_pc2))
fviz_pca_var(dades_pca, 
             select.var = list(name = rownames(var_contrib)[top_vars]),
             col.var = "contrib", # Color segons la contribució
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(title = "Contribució dels compostos més importants als components principals")



# Convertir la variable resposta a factor binari
rowData(se)$condition <- as.factor(ifelse(rowData(se)$Muscle.loss == "cachexic", 1, 0))

# Seleccionar els compostos amb major contribució basats en el PCA
top_vars <- rownames(var_contrib)[order(var_contrib[, "Dim.2"], decreasing = TRUE)[1:5]] # Canviar segons quin component principal es vol utilitzar
data_model <- as.data.frame(assay(se)[ ,top_vars]) # Seleccionem només els compostos més importants
data_model$condition <- rowData(se)$condition

# Entrenar el model de regressió logística
log_model <- glm(condition ~ ., data = data_model, family = binomial)
summary(log_model)


# Guardar les dades
save(se, file = "human_cachexia_se.Rda") # Guardar l'objecte de forma binària

# Exportar dades del SummarizedExperiment
write.csv(as.data.frame(assays(se)$counts), "human_cachexia_data.csv", row.names = FALSE)