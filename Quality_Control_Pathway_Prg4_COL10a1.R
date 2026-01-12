# --------------------------------------------------
# Spatial Transcriptomics Analysis (R version)
# Focus: Validation with PRG4 (superficial) and COL10A1 (hypertrophic)
# Universe/background = all expressed genes tested in this analysis
# --------------------------------------------------

set.seed(42)

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(glmGamPoi)
library(ggrepel)
library(mclust)
library(BayesSpace)
library(SpatialExperiment)
library(spatialLIBD)
library(GPTCelltype)
library(RColorBrewer)
library(gridExtra)
library(clusterProfiler)
library(org.Mm.eg.db) # use org.Hs.eg.db if human

# --------------------------------------------------
# 1) Load Visium Data + Expert Annotations
# --------------------------------------------------
slide1 <- Load10X_Spatial("data/slide1_WTM")

annotations <- read.csv("data/slide1_WTM/GP_features_1_Slide_1_A_WTM_Results_GP_ann.csv")
annotations <- annotations %>%
  mutate(GP = case_when(
    GP == "chondrocytes" ~ "chondrocyte",
    GP %in% c("pre-osteo","pre-osteoblasr") ~ "pre-osteoblast",
    GP == "secondary hypertophic" ~ "secondary hypertrophic",
    .default = GP
  ))

slide1$Expert_Annotation <- annotations$GP
slide1 <- slide1[, slide1$Expert_Annotation != ""] # keep annotated spots only

# Visual check
SpatialDimPlot(slide1, group.by = "Expert_Annotation", pt.size.factor = 3)

# --------------------------------------------------
# 2) Preprocessing (2000 HVGs, Resolution = 0.8)
# --------------------------------------------------
slide1 <- SCTransform(slide1, assay = "Spatial", verbose = FALSE)
slide1 <- RunPCA(slide1, assay = "SCT", verbose = FALSE)
slide1 <- FindNeighbors(slide1, dims = 1:30, k.param = 10)
slide1 <- FindClusters(slide1, resolution = 0.8)   # fixed resolution
slide1 <- RunUMAP(slide1, dims = 1:30)

SpatialDimPlot(slide1, group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Automated Clusters (Res = 0.8)")

# --------------------------------------------------
# 3) Marker Genes + Validation Zones
# --------------------------------------------------
degs <- FindAllMarkers(slide1, group.by = "Expert_Annotation")

genes_primary   <- c("Prg4")     # superficial zone
genes_secondary <- c("Col10a1")  # hypertrophic zone

SpatialFeaturePlot(slide1, features = genes_primary)   + ggtitle("PRG4 (Superficial)")
SpatialFeaturePlot(slide1, features = genes_secondary) + ggtitle("COL10A1 (Hypertrophic)")

# --------------------------------------------------
# 4) Pathway Analysis (Validation only) with explicit UNIVERSE
# --------------------------------------------------
Idents(slide1) <- "Expert_Annotation"

# DEGs (zone vs all others)
superficial_markers  <- FindMarkers(slide1, ident.1 = "superficial",  ident.2 = NULL)
hypertrophic_markers <- FindMarkers(slide1, ident.1 = "hypertrophic", ident.2 = NULL)

# ---- Define the universe/background as ALL genes expressed & tested ----
# After preprocessing, these are the genes in the Seurat object (SCT assay)
universe_genes <- unique(rownames(slide1))  # expressed/tested genes

# Foreground gene sets (padj < 0.05), intersected with the universe for safety
superficial_genes  <- intersect(rownames(superficial_markers)[superficial_markers$p_val_adj < 0.05], universe_genes)
hypertrophic_genes <- intersect(rownames(hypertrophic_markers)[hypertrophic_markers$p_val_adj < 0.05], universe_genes)

# Enrichment against the explicit universe
ego_superficial <- enrichGO(
  gene          = superficial_genes,
  universe      = universe_genes,
  OrgDb         = org.Mm.eg.db,   # swap to org.Hs.eg.db if human
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH"
)

ego_hypertrophic <- enrichGO(
  gene          = hypertrophic_genes,
  universe      = universe_genes,
  OrgDb         = org.Mm.eg.db,   # swap to org.Hs.eg.db if human
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH"
)

# Plots
dotplot(ego_superficial,  showCategory = 10) + ggtitle("Superficial Pathway Enrichment (PRG4)")
dotplot(ego_hypertrophic, showCategory = 10) + ggtitle("Hypertrophic Pathway Enrichment (COL10A1)")

# Optional: inspect counts
message("#foreground_superficial = ", length(superficial_genes),
        " | #foreground_hypertrophic = ", length(hypertrophic_genes),
        " | #universe = ", length(universe_genes))