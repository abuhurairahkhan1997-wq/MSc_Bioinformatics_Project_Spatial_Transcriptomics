# Run this in R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BayesSpace")







# 📦 Load libraries
library(SpatialExperiment)
library(BayesSpace)
library(readr)
library(S4Vectors)

# 📂 Set working directory
setwd("C:/Users/abuhu/cellpie/cellpie")

# 📥 STEP 1: Load CSV files
counts <- as.matrix(read.csv("counts_matrix.csv", header = FALSE))
coords <- as.matrix(read.csv("spatial_coords.csv", row.names = 1))
meta <- read.csv("metadata.csv", row.names = 1)
meta <- DataFrame(meta)

# ✅ Add dummy gene names BEFORE creating SpatialExperiment
gene_names <- paste0("gene_", seq_len(nrow(counts)))
rownames(counts) <- gene_names

# 🧬 STEP 2: Create SpatialExperiment object
spe <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = coords,
  rowData = DataFrame(gene_id = gene_names)  # ✅ initialize rowData here
)
colData(spe) <- meta

# 🔬 STEP 3: Select top 2000 HVGs by variance
gene_vars <- apply(counts, 1, var)
top_genes <- order(gene_vars, decreasing = TRUE)[1:2000]
rowData(spe)$is.HVG <- FALSE
rowData(spe)$is.HVG[top_genes] <- TRUE

# 🔄 STEP 4: Run BayesSpace using HVGs
spe <- spatialPreprocess(spe, n.PCs = 15)
spe <- spatialCluster(spe, q = 5, platform = "Visium", d = 15, init.method = "mclust")

# 💾 STEP 5: Save output
write.csv(colData(spe)$spatial.cluster, "bayesspace_clusters.csv")
cat("✅ BayesSpace clustering complete. Output saved as bayesspace_clusters.csv\n")

###########################################################
### Other Data set Mouselungs.

# 📦 Install and load required libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("BayesSpace", "SpatialExperiment"))

# 📦 Load libraries
library(SpatialExperiment)
library(BayesSpace)
library(readr)
library(S4Vectors)

# 📂 Set working directory to processed dataset folder
setwd("C:/Users/abuhu/cellpie/cellpie/Lungs_Mouse_Processed")

# 📥 STEP 1: Load CSV files
counts <- as.matrix(read.csv("counts_matrix.csv", header = FALSE))
coords <- as.matrix(read.csv("spatial_coords.csv", row.names = 1))
meta <- read.csv("metadata.csv", row.names = 1)
meta <- DataFrame(meta)

# ✅ Add dummy gene names BEFORE creating SpatialExperiment
gene_names <- paste0("gene_", seq_len(nrow(counts)))
rownames(counts) <- gene_names

# 🧬 STEP 2: Create SpatialExperiment object
spe <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = coords,
  rowData = DataFrame(gene_id = gene_names)
)
colData(spe) <- meta

# 🔬 STEP 3: Select top 2000 HVGs by variance
gene_vars <- apply(counts, 1, var)
top_genes <- order(gene_vars, decreasing = TRUE)[1:2000]
rowData(spe)$is.HVG <- FALSE
rowData(spe)$is.HVG[top_genes] <- TRUE

# 🔄 STEP 4: Run BayesSpace using HVGs
spe <- spatialPreprocess(spe, n.PCs = 15)
spe <- spatialCluster(spe, q = 5, platform = "Visium", d = 15, init.method = "mclust")

# 💾 STEP 5: Save clustering output
write.csv(colData(spe)$spatial.cluster, "bayesspace_clusters.csv")
cat("✅ BayesSpace clustering complete. Output saved as bayesspace_clusters.csv\n")
