library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


mc.data <- Read10X(data.dir = "/home/humebc/VTK_22/macparland/results/aggr/macparland/outs/count/filtered_feature_bc_matrix", min.cells = 3)

mc.data[1:5, 1:30]

mc <- CreateSeuratObject(counts = mc.data, project = "macparland")

mc[["percent.mt"]] <- PercentageFeatureSet(mc, pattern = "^MT-")

VlnPlot(mc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

summary(mc@meta.data)

ggplot(mc@meta.data, aes(x=nCount_RNA, y=percent.mt)) + geom_point(alpha=0.5, size=2, aes(color=nCount_RNA)) + scale_x_log10() + xlim(500, 50000) + scale_color_gradient(low="blue", high="red")

# FeatureScatter(mc, feature1 = "nCount_RNA", feature2 = "percent.mt")

# In the paper they filter by <1500 UMI and high mito genome > 0.5
mc <- subset(mc, subset = nCount_RNA > 1500 & percent.mt < 50)

VlnPlot(mc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggplot(mc@meta.data, aes(x=nCount_RNA, y=percent.mt)) + geom_point(alpha=0.5, size=2, aes(color=nCount_RNA)) + scale_x_log10() + xlim(500, 50000) + scale_color_gradient(low="blue", high="red")

mc <- NormalizeData(mc, normalization.method = "LogNormalize", scale.factor = 10000)

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
dim.heat.map = DimHeatmap(mc, dims = 1, cells = 500, balanced = T, fast=F, nfeatures=100)
ggsave("dimheatmap.png", dim.heat.map)

mc <- JackStraw(mc, num.replicate = 100)
mc <- ScoreJackStraw(mc, dims = 1:20)

JackStrawPlot(mc, dims = 1:15)

ElbowPlot(mc)

mc <- FindNeighbors(mc, dims = 1:10)
# Some further details on the parameters used for the FindCluster command can be found here:
# https://github.com/BaderLab/singleLiverCells
mc <- FindClusters(mc, resolution = 0.5, pc.use)

head(Idents(mc), 5)

mc <- RunUMAP(mc, dims = 1:10)

DimPlot(mc, reduction = "umap", pt.size=1)
ggsave("umap.png")

mc <- RunTSNE(mc, dims = 1:10)

DimPlot(mc, reduction = "tsne", pt.size=1)
ggsave("tsne.png")

saveRDS(mc, file = "mc_vtk22.rds")


# This gives us the differentially expressed genes
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
