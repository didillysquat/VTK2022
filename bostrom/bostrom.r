# All R code associated with the processing of the böstrom data

library(dplyr)
library(tidyverse)
library(GenomicFeatures)
library(RMariaDB)
library(stringr)
library(tximport)
library(DESeq2)
library("pheatmap")
library("vsn")
library(matrixStats)

txdb <- makeTxDbFromEnsembl("Homo sapiens", release=107)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

# Load the dataframe 
samples = read.csv("/home/humebc/VTK_22/bostrom/bostrom_meta.csv", header=TRUE)

# We're actually only interested in HeLa-fucci
# Filter out the U20S cells so that we only carry the HeLa cells on in the analysis
samples = samples %>% dplyr::filter(cell_type=="HeLa")
samples = samples %>% mutate(cell_type = as.factor(cell_type), cell_cycle_stage = as.factor(cell_cycle_stage), rep = as.factor(rep))

# Make a vector that contains the full paths to the abundance.h5 files
kallisto.base.dir = "/home/humebc/VTK_22/bostrom/kallisto"
files <- file.path(kallisto.base.dir, samples$dir_name, "abundance.h5")

# Verify that all the files are there
all(file.exists(files))
file.exists(files)

# Create sample names based on the meta info
samples$sample_name = stringr::str_c(samples$dir_name, "_", samples$cell_type, "_", samples$cell_cycle_stage, "_", samples$rep)
names(files) <- samples$sample_name
rownames(samples) = samples$sample_name

# Finally we can do use tximport to read in the abundance tables
txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~ cell_cycle_stage)

# Fit the model and run the DE analysis
dds = DESeq(dds)

# take a look at the results table
res = results(dds)
res_ordered = res[order(res$log2FoldChange),]

# Take a look at the other result names available 
# resultsNames(dds)
# res <- results(dds, contrast=c("cell_cycle_stage", "G2", "G1"))

# resOrdered <- res[order(res$pvalue),]
# summary(res)


# Let's perform the shrinkage so that we can see what effect it has.
resLFC <- lfcShrink(dds, coef="cell_cycle_stage_G2_vs_G1", type="apeglm")

# # Let's compare the shrinkage results
# plotMA(res, ylim=c(-2,2))
# plotMA(resLFC, ylim=c(-2,2))

# Two forms of normalization to do the plotting with
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)

# Finally we can get on to the heat map

# ENSG00000226964 is -8.7 logfold change
# ENSG00000237599 is 21.7 logfold change


# I don't know how they have ordered the cells in the figure they made but fist I'd like to try
# going in order of the log2 fold change S compared to G1
res = results(dds)
res_ordered = res[order(res$log2FoldChange, decreasing=TRUE),]
head(res)

df <- as.data.frame(colData(dds)[,c("cell_cycle_stage")])
colnames(df) = c("cell_cycle_stage")
rownames(df) = rownames(samples)
select <- order(res$log2FoldChange, decreasing=FALSE)[1:10000]
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, scale="row")


# Let's say we're only interested in those genes that were differentially expressed between one of the
# 3 comparisons
# We can generate a list of these genes
s1_g1_sig_res = as.data.frame(results(dds, contrast=c("cell_cycle_stage", "S", "G1"))) %>% dplyr::filter(padj <= 0.01)
dim(s1_g1_sig_res)
s1_g2_sig_res = as.data.frame(results(dds, contrast=c("cell_cycle_stage", "S", "G2"))) %>% dplyr::filter(padj <= 0.01)
dim(s1_g2_sig_res)
g1_g2_sig_res = as.data.frame(results(dds, contrast=c("cell_cycle_stage", "G1", "G2"))) %>% dplyr::filter(padj <= 0.01)
dim(g1_g2_sig_res)

sig_genes_redundant = c(rownames(s1_g1_sig_res), rownames(s1_g2_sig_res), rownames(g1_g2_sig_res))
length(sig_genes_redundant) # 4129
sig_genes_unique = unique(sig_genes_redundant)
length(sig_genes_unique) # 2677

# Now we have a set of genes that are differentiatlly expressed in one or more of the comparisons

# Next we need to make a df where we have the average count by cell_cycle_stage
dds_counts = as.data.frame(assay(dds))
dds_counts_w_av = dds_counts %>% dplyr::mutate(
    s_av=rowMeans(select(., c(SRR6150372_HeLa_S_1, SRR6150373_HeLa_S_2, SRR6150374_HeLa_S_3))),
    g1_av=rowMeans(select(., c(SRR6150369_HeLa_G1_1, SRR6150370_HeLa_G1_2, SRR6150371_HeLa_G1_3))),
    g2_av=rowMeans(select(., c(SRR6150375.combined_HeLa_G2_1, SRR6150378.combined_HeLa_G2_2, SRR6150381.combined_HeLa_G2_3)))
    )
# Finally we want to filter down to just the average columns
dds_counts_only_av = dds_counts_w_av %>% select(c(s_av, g1_av, g2_av))
dds_counts_only_av_sig = dds_counts_only_av[sig_genes_unique,]

# Then given that there was complete overlap, we can plot the genes in order of abundance
# but first we will want to normalize the data by row so that we can pre-order by post
# normalized abundance.
head(dds_counts_only_av_sig)
dds_counts_only_av_sig_s1_ordered = dds_counts_only_av_sig[order(dds_counts_only_av_sig$s_av, decreasing = T),]
pheatmap(dds_counts_only_av_sig_s1_ordered[1:10,],cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE)


RowSD <- function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}
# We need to standardize the data by row and then order and then plot
dds_counts_only_av_sig_stand = dds_counts_only_av_sig %>% 
  mutate(mean = (s_av + g1_av + g2_av)/3, stdev = RowSD(cbind(s_av, g1_av, g2_av)))

dds_counts_standardized = dds_counts_only_av_sig_stand %>% mutate(
    s_std = (s_av-mean)/stdev,
    g1_std = (g1_av-mean)/stdev,
    g2_std = (g2_av-mean)/stdev
    ) %>% dplyr::select(s_std, g1_std, g2_std) %>% arrange(desc(s_std))

pheatmap(dds_counts_standardized, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE)

# At this point I officially gave up! They just don't seem to agree.