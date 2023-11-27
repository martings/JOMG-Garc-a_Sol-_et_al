# Load necessary libraries
library(Seurat)
library(scales)

# Read in the Seurat object from an RDS file
#pbmc <- readRDS("path/to/your/file/")
pbmc <- readRDS("/media/megas/0123-4567/single_cell_reference/Integracion_solo_post_natal_Back_Up_8_9.rds")

# Plot UMAP visualization of the Seurat object
DimPlot(pbmc, reduction = "umap")

# Feature plot for the "Procr" gene with cutoffs
FeaturePlot(pbmc, features = "Procr", min.cutoff = "p10", max.cutoff = "p90", slot = "data")

# Read in a subset of the data from another RDS file
#pbmc_sub <- readRDS("path/to/your/file/Subset_Luminales_Alveolares_pbmc_seurat.rds")
pbmc_sub <- readRDS("/media/megas/0123-4567/single_cell_reference/Subset_Luminales_Alveolares_pbmc_seurat.rds")

# UMAP visualization of the Seurat subset by cell type
DimPlot(pbmc, reduction = "umap", group.by = "Cell_Type")

# UMAP visualization of the subset by cell type
DimPlot(pbmc_sub, reduction = "umap", group.by = "Cell_Type")

# Show color palette for the 'Cell_Type' variable
show_col(hue_pal()(12))

# UMAP visualization with custom colors for different cell types
DimPlot(object = pbmc_sub, reduction = 'umap', group.by = 'Cell_Type',
        cols = c('Luminal_Progenitor' = '#F8766D', 'Luminal_AlvProg' = '#00BFC4', 'Luminal_AlvSec' = '#00B4F0'))

# UMAP visualization focusing on a specific cell type ('Luminal_AlvSec') with a custom color
DimPlot(object = pbmc_sub, reduction = 'umap', group.by = 'Cell_Type',
        cols = c('Luminal_AlvSec' = '#00BFC4'))

# Blended feature plots for selected genes
FeaturePlot(pbmc, features = c("Zfp36", "Aldh1a3"), blend = TRUE, min.cutoff = "p10", max.cutoff = "p90", pt.size = 1)

FeaturePlot(pbmc, features = c("Zfp36", "Procr"), blend = TRUE, min.cutoff = "p10", max.cutoff = "p90", pt.size = 1)
