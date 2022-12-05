# All R code associated with the processing of the böstrom data

library(dplyr)
library(tidyverse)
library(GenomicFeatures)
library(RMariaDB)
library(stringr)
library(tximport)
library(DESeq2)

txdb <- makeTxDbFromEnsembl("Homo sapiens", release=107)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Load the dataframe 
samples = read.csv("/home/humebc/VTK_22/bostrom/bostrom_meta.csv", header=TRUE)

# We're actually only interested in HeLa-fucci
# Filter out the U20S cells so that we only carry the HeLa cells on in the analysis
samples = samples %>% dplyr::filter(cell_type=="HeLa")
samples = samples %>% mutate(cell_type = as.factor(cell_type),cell_cycle_stage = as.factor(cell_cycle_stage), rep = as.factor(rep))

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



