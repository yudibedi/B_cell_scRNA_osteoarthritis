# Install the required packages----

install.packages("BiocManager")
install.packages('Seurat')
install.packages("devtools")
remotes::install_github("mojaveazure/seurat-disk")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scran")
BiocManager::install("scater")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")
install.packages("vegan")
install.packages("pheatmap")
install.packages("tidyverse")
install.packages("enrichR")
install.packages("scCustomize")

sessionInfo()

library(Seurat)
library(SeuratDisk)
library(SingleR)
library(celldex)
library(dplyr)
library(Matrix)
library(scran)
library(scater)
library(vegan)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(enrichR)
library(ggpubr)
library(scCustomize)

rm(list = ls())
# setting location paths
data_path <- "~/data/"
fig_path <- "~/figures/"

## Convenience functions ----
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}

# Reading in data ----

## reading count matrix generated using Parse pipeline
mat_path <- "~/data/DGE_filtered"
mat <- ReadParseBio(mat_path)

## Check to see if empty gene names are present, add name if so.

table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data

cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), 
                      row.names = 1)

# Create object

pbmc <- CreateSeuratObject(mat, 
                           min_genes = 100, 
                           min_cells = 100, 
                           names.feild = 0, 
                           meta.data = cell_meta)

# Setting our initial cell class to a single type, this will change after clustering. 

pbmc@meta.data$orig.ident <- factor(rep("pbmc", 
                                        nrow(pbmc@meta.data)))
Idents(pbmc) <- pbmc@meta.data$orig.ident

# save the object up to this point so it can be loaded quickly

SaveObject(pbmc, "seurat_obj_before_QC")
pbmc <- ReadObject("seurat_obj_before_QC")


## Cell quality control ----
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
plot <- VlnPlot(pbmc, pt.size = 0.10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SaveFigure(plot, "vln_QC", width = 12, height = 6)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
SaveFigure((plot1 + plot2),"scatter_QC", width = 12, height = 6, res = 200)

# Perform the filtering for outliers

pbmc <- subset(pbmc, 
               subset = nFeature_RNA < 5000 & 
                 nCount_RNA < 20000 & 
                 percent.mt < 15)

## Normalizing, scaling and PCA reduction of the data ----

pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

## Identification of highly variable features

pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

### Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(pbmc), 10)

### Plot variable features with and without labels

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
SaveFigure((plot1 + plot2), "var_features", width = 12, height = 6)

## Scaling the data to apply linear transformation

pbmc <- ScaleData(pbmc)

## perform linear dimension reduction

pbmc <- RunPCA(pbmc)

## save the object up to this point so it can be loaded quickly

SaveObject(pbmc, "seurat_obj_after_PCA")

pbmc <- ReadObject("seurat_obj_after_PCA")

# Examine and visualize PCA results a few different ways - VizDimloading, DimPlot and DimHeatmap

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

plot <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
SaveFigure(plot, "viz_PCA_loadings", width = 10, height = 8)

plot <- DimPlot(pbmc, reduction = "pca", group.by = "orig.ident")
SaveFigure(plot, "pc1_2_scatter", width = 8, height = 6)

# DimHeatmap - Image doesn't save as png unless fast = FALSE

plot <- DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
SaveFigure(plot, "dim_heatmap1_15", width = 12, height = 18)

# Determine the dimensionality of the dataset
# NOTE: Jackstraw plots can take a long time for big datasets. 
# More approximate techniques such as those implemented in ElbowPlot() 
# can be used to reduce computation time

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
plot <- JackStrawPlot(pbmc, dims = 1:20)
SaveFigure(plot, "PC_jackstraw_plot", width = 10, height = 8)

# alternate method is to generate elbow plot
plot <- ElbowPlot(pbmc,ndims = 50)
SaveFigure(plot, "PC_elbow_plot", width = 8, height = 10)

# Clustering  ----

# cluster the cells - these values were chosen to make 15 dimensions of reduction
# and resolution of 0.6
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.60)

# reordering clusters according to their similarity - 
# done for ease in merging cell clusters by reassigning each cluster position in phylogenic tree

install.packages("ape")
library(ape)

pbmc <- BuildClusterTree(pbmc, 
                         reorder = TRUE, 
                         reorder.numeric = TRUE)

## Run UMAP non-linear dimensional reduction ----

pbmc <- RunUMAP(pbmc, dims = 1:15)
plot <- DimPlot(pbmc, reduction = "umap")
SaveFigure(plot, "umap_louvain_res_06", width = 9, height = 8)

# Sanity check
pbmc

# Saving the object at this point so it can be easily loaded back in without 
# computationally intensive steps or to share with collaborators

SaveObject(pbmc, "seurat_obj_clustered")
pbmc <- ReadObject("seurat_obj_clustered")

## Finding cluster markers based on diff gene expression ----

pbmc_markers <- FindAllMarkers(pbmc, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

### Visualizing top n genes per cluster in a single figure ----
# DoHeatmap is a good tool for few thousand cell but loses resolution as number of cells increase
# DotPlot addresses this by displaying average expression for all cells in a cluster
top5 <- pbmc_markers %>% 
  group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot <- unique(top5$gene)
plot <- DotPlot(pbmc, 
                features = to_plot, 
                group.by = "tree.ident") + coord_flip()
SaveFigure(plot, "dplot_top5", width = 9, height = 20)

# visualizing gene marker expression  by VlnPlot - B cell markers
# this should be a bulk of our cell type as we enriched PBMCs for B cells
plot <- VlnPlot(pbmc, 
                features = c("PAX5", "MS4A1", "CD79A"), 
                group.by = "tree.ident")
SaveFigure(plot, "vln_B_markers", width = 16, height = 8)

# you can plot raw counts as well - NK cell markers
plot <- VlnPlot(pbmc, 
                features = c("NKG7", "GNLY"), 
                slot = "counts", 
                log = TRUE)
SaveFigure(plot, "vln_NK_cell_markers", width = 16, height = 8)

# Cluster classsification ----

#  Violin plot with known PBMC markers from Seurat tutorial
plot <- VlnPlot(pbmc, 
                features = c("MS4A1", "GNLY", "CD3E", 
                             "CD14", "FCER1A", "FCGR3A", 
                             "LYZ", "PPBP","CD8A"), 
                group.by = "tree.ident")
SaveFigure(plot, "vln_exp_known_PBMC_markers", width = 48, height = 32)

#  Violin plot with some B cell markers from papers
plot <- VlnPlot(pbmc, features = c("PAX5", "MS4A1", "CD79A", "IGHM", "IGHD", "TBX21", "LTB", "ITGAX", "IL32"), group.by = "tree.ident")
SaveFigure(plot, "vln_exp_new_markers", width = 16, height = 8)

#  Violin plot with B cell markers to distinguish swithced memory B from exhausted B
plot <- VlnPlot(pbmc, features = c("MS4A1", "CD79A", "CD27", "TIGIT", "PDCD1", "CD160", "PDCD1LG2", "CEACAM21", "BACH2"), group.by = "tree.ident", ncol = 4) 
SaveFigure(plot, "vln_exp_new_B_markers", width = 16, height = 8)


## Automated classification using  singleR ----
# the following is from BIG Bioinformatics scRNA workshop 
# https://www.bigbioinformatics.org/intro-to-scrnaseq

# Predict using the MonacoImmuneData of bulk RNA markers of PBMCs

mimd.se <- MonacoImmuneData()

# What is the structure of mimd.se?
mimd.se
View(as.data.frame(mimd.se@colData))

# Predict using the Broad labels
sceP <- as.SingleCellExperiment(pbmc)
pred.mimd <- SingleR(test = sceP, 
                     ref = mimd.se, 
                     labels = mimd.se$label.main)
plot <- plotScoreHeatmap(pred.mimd,
                         max.labels = 8, 
                         clusters = pbmc@active.ident, 
                         order.by = "clusters", 
                         show_colnames = FALSE)
SaveFigure(plot, "broad_label_markers", width = 10, height = 6)

# Predict using the Fine-grained labels
Idents(pbmc) <- pbmc$tree.ident
pred.mimd.fine <- SingleR(test = sceP, 
                          ref = mimd.se, 
                          labels = mimd.se$label.fine)
plot <- plotScoreHeatmap(pred.mimd.fine,
                         max.labels = 14,
                         clusters = pbmc@active.ident,
                         order.by = "clusters",
                         show_colnames = FALSE)
SaveFigure(plot, "fine_label_markers", width = 10, height = 15)


# Add the predicted labels and compare to clusters
# You are looking for clusters which have unambiguous labels

pred.mimd$pruned.labels
pred.mimd.fine$pruned.labels

pbmc$predicted_broad_id <- pred.mimd$pruned.labels
pbmc$predicted_fine_id <- pred.mimd.fine$pruned.labels

table(pbmc$predicted_fine_id, pbmc@active.ident)
table(pbmc$predicted_broad_id, pbmc@active.ident)

## Selecting final IDs ----
new_ids <- c(
  "1" = "Exhausted B cells",
  "2" = "Classical monocytes",
  "3" = "Intermediate monocytes",
  "4" = "Progenitor cells",
  "5" = "Naive B cells",
  "6" = "Naive B cells",
  "7" = "Switched memory B cells",
  "8" = "Naive B cells",
  "9" = "Non-switched memory B cells",
  "10" = "Natural killer cells",
  "11" = "T cells",
  "12" = "Dendritic cells",
  "13" = "Naive B cells",
  "14" = "Naive B cells"
)

pbmc <- RenameIdents(pbmc, new_ids)
pbmc$cell_type <- Idents(pbmc)

plot <- DimPlot_scCustom(pbmc, 
                         label.size = 5, 
                         label = TRUE,
                         repel = TRUE, 
                         pt.size = 0.1) + 
  NoLegend()

# side-by-side view of seurat clusters, samples and final IDs
Idents(pbmc) <- pbmc$cell_type

plot <- DimPlot(pbmc, group.by = "tree.ident", 
                label.size = 5, 
                label = TRUE,
                repel = TRUE) +  NoLegend() +
  DimPlot(pbmc, 
          group.by = "sample") +
  DimPlot(pbmc, 
          label.size = 5, 
          label = TRUE,
          repel = TRUE ) + NoLegend()

SaveFigure(plot, "pbmc_dimplot_comparison", width = 20, height = 6)


# Saving the object at this point so it can be easily loaded back in without 
# computationally intensive steps or to share with collaborators
SaveObject(pbmc, "seurat_obj_clustered_new_IDs")
