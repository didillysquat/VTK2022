library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Here we read in the filtered features matrix that was output from the Cell Ranger aggr analysis
mc.data <- Read10X(data.dir = "/home/humebc/VTK_22/macparland/results/aggr/macparland/outs/count/filtered_feature_bc_matrix", min.cells = 3)

# Have a look at the count matrix
mc.data[1:5, 1:30]

# Create a Seurat Object
mc <- CreateSeuratObject(counts = mc.data, project = "macparland")

# Create the percent mitochondrial reads mapped attribute
mc[["percent.mt"]] <- PercentageFeatureSet(mc, pattern = "^MT-")

# Create a violin plot to look at each of the features
VlnPlot(mc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Summarise the meta data
summary(mc@meta.data)

# Try to create the first figure with the built in method
FeatureScatter(mc, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Have a little go at making it similar
ggplot(mc@meta.data, aes(x=nCount_RNA, y=percent.mt)) + geom_point(alpha=0.5, size=2, aes(color=nCount_RNA)) + scale_x_log10() + xlim(500, 50000) + scale_color_gradient(low="blue", high="red")

# Why don't they agree?

# In the paper they filter by <1500 UMI and high mito genome > 0.5
mc <- subset(mc, subset = nCount_RNA > 1500 & percent.mt < 50)

# Verify with the violin plot that the filtering has taken effect
VlnPlot(mc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Plot up the ggplot version of the figure after filtering
ggplot(mc@meta.data, aes(x=nCount_RNA, y=percent.mt)) + geom_point(alpha=0.5, size=2, aes(color=nCount_RNA)) + scale_x_log10() + xlim(500, 50000) + scale_color_gradient(low="blue", high="red")

# They still don't agree, right?

# Standard Seurat normalization acts by dividing by the total number of counts,
# multiplies by a scale factor (default 10000) the log normalizes
mc <- NormalizeData(mc, normalization.method = "LogNormalize", scale.factor = 10000)

# Generally researchers work with high variable features for the continuing analysis
mc <- FindVariableFeatures(mc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate

# Also referred to as standardization
all.genes <- rownames(mc)
mc <- ScaleData(mc, features = all.genes)

# Perform PCA
mc <- RunPCA(mc, features = VariableFeatures(object = mc))

VizDimLoadings(mc, dims = 1:2, reduction = "pca")
ggsave("visdimloadings.png")

DimPlot(mc, reduction = "pca")
ggsave("dimplot.png")

# NB fast=F must be set in order to return a ggplot object.
DimHeatmap(mc, dims = 1, cells = 500, balanced = T, fast=F, nfeatures=100)
DimHeatmap(mc, dims = 1:15, cells = 500, balanced = TRUE, fast=F)
ggsave("dimheatmap.png")

mc <- JackStraw(mc, num.replicate = 100)
mc <- ScoreJackStraw(mc, dims = 1:20)

JackStrawPlot(mc, dims = 1:15)

# The authors used 29 PCs
# Take a look and see if that was an appropriate call
ElbowPlot(mc, ndims=35)


# Create the graph to perform the smart clustering on.
mc <- FindNeighbors(mc, dims = 1:10)

# Some further details on the parameters used for the FindCluster command can be found here:
# https://github.com/BaderLab/singleLiverCells
mc <- FindClusters(mc, resolution = 0.5, pc.use=1:10)

head(Idents(mc), 5)

# Make a UMAP
mc <- RunUMAP(mc, dims = 1:10)

DimPlot(mc, reduction = "umap", pt.size=1)
ggsave("umap.png")

# Make a tSNE plot
mc <- RunTSNE(mc, dims = 1:10)

DimPlot(mc, reduction = "tsne", pt.size=1)
ggsave("tsne.png")

saveRDS(mc, file = "mc_vtk22.rds")


# This gives us the differentially expressed genes for each of the clusters
# That can then be used as marker genes to identify the clusters.
mc.markers <- FindAllMarkers(mc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


mc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(mc, features = top10$gene) + NoLegend()
ggsave("cluster_diff_expression.png")

# From here you would have to then work with their curated list of markers
# This is rather inaccessible and a lot of work so we will call it quits at this point.
# Nice work!
