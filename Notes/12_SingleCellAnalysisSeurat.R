library(dplyr)
library(Seurat)
library(patchwork)

### For more details, please check: https://satijalab.org/seurat/articles/pbmc3k_tutorial



### Accessing the data.

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/storage/ice-shared/biol6150/Data/SingleCell/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

#2,700 single cells that were sequenced on the Illumina NextSeq 500
dim(pbmc.data)
head(pbmc.data)

#Sparse vs Dense matrix.
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

sparse.size <- object.size(pbmc.data)
sparse.size


### Quality control.
# 1. Low-quality cells or empty droplets will often have very few genes
# 2. Cell doublets or multiplets may exhibit an aberrantly high gene count
# 3. Low-quality / dying cells often exhibit extensive mitochondrial contamination

# They calculate mitochondrial QC metrics with the PercentageFeatureSet() function, 
#which calculates the percentage of counts originating from a set of features

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

## QC Approach:

# They visualize QC metrics, and use those to filter cells.
# 1. They filter cells that have unique feature counts over 2,500 or less than 200
# 2. They filter cells that have >5% mitochondrial counts

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Pruning
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

### Normalization
# method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

### Feature selection
# highly expressed in some cells, and lowly expressed in others
# models on mean-variance relationship inherent in single-cell data
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


### Scaling the data. (linear transformation)
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


### PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#Check if cells cluster in PCA?
DimPlot(pbmc, reduction = "pca") + NoLegend()

# Allows for easy exploration of the primary sources of heterogeneity in a dataset, 
# and can be useful when trying to decide which PCs to include for further downstream analyses. 
# Both cells and features are ordered according to their PCA scores
DimHeatmap(pbmc, dims = 1:3, cells = 500, balanced = TRUE)

### Determine dimensionality of the dataset.
# Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one 
ElbowPlot(pbmc)

### Clustering and visualization.
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
# Recommended not to draw biological conclusions from visualization techniques.

# Clustering using a graph based approach: KNN-graph.
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

### Markers specific to each cluster.
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the significant ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)



































