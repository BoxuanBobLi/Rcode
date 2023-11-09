library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)

Psoriasis <- readRDS("/Users/liboxuan/Desktop/学习/biomathU/data/Psoriasis.rds")

Psoriasis[["percent.mt"]] <- PercentageFeatureSet(Psoriasis, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Psoriasis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Psoriasis, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Psoriasis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Psoriasis <- subset(Psoriasis, subset = nFeature_RNA > 200 & percent.mt < 20)
#noamalization
Psoriasis <- NormalizeData(Psoriasis, normalization.method = "LogNormalize", scale.factor = 10000)
Psoriasis <- NormalizeData(Psoriasis)

#Identification of highly variable features (feature selection)
Psoriasis <- FindVariableFeatures(Psoriasis, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Psoriasis), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Psoriasis)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#Scaling the data
all.genes <- rownames(Psoriasis)
Psoriasis <- ScaleData(Psoriasis, features = all.genes)

#Perform linear dimensional reduction
Psoriasis <- RunPCA(Psoriasis, features = VariableFeatures(object = Psoriasis))
# Examine and visualize PCA results a few different ways
print(Psoriasis[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Psoriasis, dims = 1:2, reduction = "pca")
DimPlot(Psoriasis, reduction = "pca")
DimHeatmap(Psoriasis, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Psoriasis, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Psoriasis <- JackStraw(Psoriasis, num.replicate = 20)
Psoriasis <- ScoreJackStraw(Psoriasis, dims = 1:20)
JackStrawPlot(Psoriasis, dims = 1:20)
ElbowPlot(Psoriasis)

#Cluster the cells
Psoriasis <- FindNeighbors(Psoriasis, dims = 1:20)
Psoriasis <- FindClusters(Psoriasis, resolution = 2)

# Look at cluster IDs of the first 5 cells
head(Idents(Psoriasis), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Psoriasis <- RunUMAP(Psoriasis, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Psoriasis, reduction = "umap")

Idents(Psoriasis) <- Psoriasis$full_clustering
cellchatPsoriasis <- createCellChat(object = Psoriasis)
cellchatPsoriasis@idents

Idents(Healthy) <- Healthy$full_clustering
cellchatHealthy <- createCellChat(object = Healthy)
cellchatHealthy@idents

