# 📦 Load Required Libraries
library(BayesSpace)
library(SpatialExperiment)
library(readr)

# 📁 Load Input Data
counts <- as.matrix(read_csv("bayesspace_input/counts_matrix.csv", col_names = FALSE))
coords <- read_csv("bayesspace_input/spatial_coords.csv")
meta <- read_csv("bayesspace_input/metadata.csv")

# 🧱 Create SpatialExperiment Object
spe <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = as.matrix(coords[, 2:3]),
  colData = meta
)
rownames(spe) <- paste0("Gene", seq_len(nrow(counts)))
colnames(spe) <- meta$Barcode

# 🔄 Preprocessing
spe <- spatialPreprocess(spe, platform = "Visium", log.normalize = TRUE)

# ⚙️ Run KMeans Initialization
set.seed(42)
spe_kmeans <- spatialCluster(spe, q = 5, d = 15, init.method = "kmeans", model = "t", nrep = 50000)
write.csv(colData(spe_kmeans)$spatial.cluster, "bayesspace_clusters_kmeans.csv", row.names = FALSE)

# ⚙️ Run Mclust Initialization
set.seed(42)
spe_mclust <- spatialCluster(spe, q = 5, d = 15, init.method = "mclust", model = "t", nrep = 50000)
write.csv(colData(spe_mclust)$spatial.cluster, "bayesspace_clusters_mclust.csv", row.names = FALSE)

# ✅ Done: Both clustering results are saved in your working directory

########################################################################################
# 📦 Load Required Libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install(c("BayesSpace", "SpatialExperiment"))

# MouseLungsdataset

library(BayesSpace)
library(SpatialExperiment)
library(readr)
library(S4Vectors)

# 📁 Load Input Data
setwd("C:/Users/abuhu/cellpie/cellpie/Lungs_Mouse_Processed")

counts <- as.matrix(read_csv("counts_matrix.csv", col_names = FALSE))
coords <- read_csv("spatial_coords.csv")
meta <- read_csv("metadata.csv")

# ✅ Assign proper gene and spot names
rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
colnames(counts) <- meta$Barcode
meta <- DataFrame(meta)

# 🧱 Create SpatialExperiment Object
spe <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = as.matrix(coords[, c("row", "col")]),
  colData = meta,
  rowData = DataFrame(gene_id = rownames(counts))
)

# 🔄 Preprocessing
spe <- spatialPreprocess(spe, platform = "Visium", log.normalize = TRUE)

# ⚙️ Run KMeans Initialization
set.seed(42)
spe_kmeans <- spatialCluster(
  spe, q = 5, d = 15, init.method = "kmeans",
  model = "t", nrep = 50000
)
write.csv(colData(spe_kmeans)$spatial.cluster, "bayesspace_clusters_kmeans.csv", row.names = FALSE)

# ⚙️ Run Mclust Initialization
set.seed(42)
spe_mclust <- spatialCluster(
  spe, q = 5, d = 15, init.method = "mclust",
  model = "t", nrep = 50000
)
write.csv(colData(spe_mclust)$spatial.cluster, "bayesspace_clusters_mclust.csv", row.names = FALSE)

cat("✅ BayesSpace clustering complete. Outputs saved to CSV files.\n")
