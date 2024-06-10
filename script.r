library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

dir <- "~/Bionformatics_Projects/HumanPBMC_Analysis"
setwd(dir)

# Load the expression matrix
adj.matrix <- Read10X(data.dir = "filtered_feature_bc_matrix/")

# Create seurat object
srat <- CreateSeuratObject(adj.matrix, project = "HumanPBMC")
srat

# Erase adj.matrix from memory to save RAM
adj.matrix <- NULL
str(srat)

# Obtaining Meta.data
meta <- srat@meta.data
head(meta)

#view
head(meta)

# Generate summary statistics
summary(meta$nCount_RNA)

# Identify mitochondrial genes
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")

# Make Violin Plots of the selected metadata features
vln_plot <- VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA" ,"percent.mt"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size = 10))
# Save the plot
ggsave("Plots/vln_plot.png", plot = vln_plot, width = 10, height = 5)

# Generate a scatter plot to visualize relationship between nCount_RNA and percent.mt
scatter_plot <- FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
# Save the plot
ggsave("Plots/scatter_plot.png", plot = scatter_plot, width = 10, height = 5)

# Generate a scatter plot to visualize relationship between nCount_RNA and nFeature_RNA
scatter_plot_nFeature_RNA <- FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Save the plot
ggsave("Plots/scatter_plot_nFeature_RNA.png", plot = scatter_plot_nFeature_RNA, width = 10, height = 5)

table(srat[['QC']])
### NORMALIZATION AND DIMENSIONALITY REDUCTION
# Normalize the data to account for sequencing depth
srat <- NormalizeData(srat)

# The most variable features (genes)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)

# To identify the 50 most highly variable genes:
top10 <- head(VariableFeatures(srat), 10)
top10

# Plot Variable features
varialble_plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = varialble_plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
# If saving to a file, specify dimensions
ggsave("Plots/variable_plot.png", plot = varialble_plot1, width = 10, height = 8)
dev.off()

# Centering and Scaling data matrix
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)

# Run Principal Component Analysis
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
# Visualization of the loadings of the PCA
pca_plot <- VizDimLoadings(srat, dims = 1:9, reduction = "pca") & theme(axis.text = element_text(size = 5), axis.title = element_text(size = 8, face = "bold"))
ggsave("Plots/pca_loadings_plot.png", plot = pca_plot, width = 10, height = 8)
dev.off()

# Heatmap of each Principal Component
pca_heatmap <- DimHeatmap(srat, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
ggsave("Plots/pca_heatmap.png", plot = pca_heatmap, width = 15, height = 12)

# Visualise dimensionality reduction
dim_plot <- DimPlot(srat, reduction = "pca")
ggsave("Plots/dim_plot.png", plot = dim_plot, width = 15, height = 12)

# Generate an elbow plot to determine the significant principal components
elbow_plot <- ElbowPlot(srat)
ggsave("Plots/elbow_plot.png", plot = elbow_plot, width = 15, height = 12)

### CLUSTERING ###
# Find Neighbours (Construct the SNN Graph)
srat <- FindNeighbors(srat, dims = 1:11)

# Cluster Identification: Find clusters
srat <- FindClusters(srat, resolution = 0.5)

# Visualise Clusters
srat <- RunUMAP(srat, dims = 1:11, verbose = F)

# Cluster sizes
table(srat@meta.data$seurat_clusters)

cluster_plot <- DimPlot(srat, label.size = 4, repel = T, label = T)
ggsave("Plots/cluster_plot.png", plot = cluster_plot, width = 15, height = 12)

# visualize the expression of four genes:  visualizing the spatial distribution of gene expression in reduced dimensionality space
expression_plot <- FeaturePlot(srat, features = c("HBB", "C1QA", "C1QB", "IFI27"))
ggsave("Plots/expression_plot.png", plot = expression_plot, width = 15, height = 12)

# Percentage of RNA expression
RNAfeatures_plot <- FeaturePlot(srat, features = "nFeature_RNA") & theme(plot.title = element_text(size = 10))
ggsave("Plots/RNAfeatures_plot.png", plot = RNAfeatures_plot, width = 15, height = 12)
# Percentage of mitochondrial gene expression
mt_features <- VlnPlot(srat, features = "percent.mt") & theme(plot.title = element_text(size = 10))
ggsave("Plots/mt_features_plot.png", plot = mt_features, width = 15, height = 12)

VlnPlot(srat, features = "percent.mt") & theme(plot.title = element_text(size = 10))

# Read in the expression matrix

# Load the updated cell cycle gene lists provided by Seurat
cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Assign Cell-Cycle Scores
srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes)

# View cell cycle scores and phase assignments
head(srat@meta.data)
table(srat@meta.data$Phase)


FeaturePlot(srat,features = "percent.mt",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size = 10))

VlnPlot(srat,features = "percent.mt") & theme(plot.title = element_text(size=10))


VlnPlot(srat,features = c("nCount_RNA","nFeature_RNA")) & theme(plot.title = element_text(size=10))


FeaturePlot(srat, features = c("S.Score", "G2M.Score"), label.size = 4, repel = T, label = T) & theme(plot.title = element_text(size=10))


VlnPlot(srat, features = c("S.Score", "G2M.Score")) & theme(plot.title = element_text(size=10))






dev.off()
