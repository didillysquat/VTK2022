# Advanced Technologies for the Life Sciences: Bioinformatics Practical
The markdown document will cover all aspects of the bioinformatics
practical related to the 2022 VTK Advanced Technologies for the Life Sciences.

# Introduction
In this practical we will process 2 RNA sequencing (RNA-seq) datasets
to generate a selection of the analytical results presented in the two
papers associated with the datasets. One of the datasets is a bulk-RNA
dataset while the other is a single cell RNA (scRNA-seq) dataset.

[What is bulk RNA-seq?](https://www.scdiscoveries.com/support/what-is-bulk-rna-sequencing/)

[What is single cell RNA-seq?](https://en.wikipedia.org/wiki/Single-cell_transcriptomics)

## The papers
These are the two papers we will be working with:

- Bulk RNA-seq: [Böstrom et al 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0188772). Comparative cell cycle
transcriptomics reveals synchronization of developmental
transcription factor networks in cancer cells. PLOS ONE.

- scRNA-seq: [MacParland et al 2018](https://www.nature.com/articles/s41467-018-06318-7). Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations. Nature Communications.

Both papers have been uploaded to the ILLIAS system in Publications/SequAna

It would be a good idea to make yourself familiar with the papers before starting the practical work.

## Location
The practical will take place in room M739 from 09.00-16.30 on the 6-8th of December. It is currently unclear whether we will have access to the room on the 9th of December as well.

Ben Hume, the current SequAna bioinformatician will
be running the practical and there to assist you.

## Computing setup
You will need access to a computer to complete this practical.
While M739 does contain computers that you could use to login to the
SequAna computation server (where much of the work will occur)
this is a very inconvenient way to work as your home directory is
deleted after every session and it is not possible to install
applications that use Graphical User Interfaces (GUI).

As such it is strongly recommended that you bring a laptop with
you to complete the practical. There will be 1 spare laptop
available for use during the practical session in M739, but this
laptop cannot leave the room as it belongs to SequAna.

# Objectives

The main objective of this practical is to give you an introduction to the tools that are used by computational biologist/bioinformaticians to generate meaningful results from sequencing data.

The objective of this course is not for you to become proficient or masterful of the techniques we will be covering (we have only 3 or 4 days!). Any proficiency gained, however, 
will likely be extremely valuable to you in your career as a research scientist.

To achieve this objective we will work with the sequencing data archived as part of the above mentioned studies to recapitulate several of their key findings.

In doing so we will cover many broad informatic/bioinformatic techniques not limited to:

- Working on the command line interface (CLI)
- Using Conda environments to install programs and packages
- Working with Docker images in Singularity
- Working with core bioinformatic tools to perform:
    - access of archived sequencing data
    - sequencing pre-processing and quality control
    - sequence analysis
- Workflow management with Nextflow
- R scripting to manipulate, analyze and visualize data

I will provide resources for all topics we cover and you are encouraged to look at these
resources if you wish to further your knowledge of the topic.

# Structure of the practical
The practical will be divided up by days (1-4).

Each day we will work towards our end goal of recapitulating the results of our chosen studies.

One of the most important skills in computation biology / informatics is the effective
sourcing of reference material. I.e. good googleing!

As such throughout the 4 days, while you will be given a structure to follow,
you will also be asked to work out how to do certain tasks on your own.
But don't worry, the SequAna bioinformatician will be there to help you if you get stuck.

If there are any requirements for the day's exercises these will be listed at the beginning of each day's section.

# DAY 1: Installing programs and fetching data

## Requirements
You should already be familiar with working in bash on the command line.
See [here](https://rnabio.org/module-00-setup/0000/08/01/Unix/
) for an excellent resource.

If you don't yet have access to a unix-based operating system you can
use [this](https://cocalc.com/projects?anonymous=terminal).

## Part 1: Connecting to SequAna's computational server: sequana
We will perform much of our analyses on SequAna's
computational server that is creatively named 'sequana'.

The server runs a Linux distribution: Ubuntu.

[What is an operating system](https://en.wikipedia.org/wiki/Operating_system)

[What is Unix](https://en.wikipedia.org/wiki/Unix)

[What is Linux?](https://www.linux.com/what-is-linux/)

The vast majority of bioinformatic programs run in Linux,
although many can also be run in Windows or on Mac OS X.

Max OS X is Unix-based and therefore has many similarities to Linux distributions.
The terminal app on Mac OS X offers the user a CLI very similar to that of Linux distributions.

To connect to sequana we will use ssh.

[What is ssh?](https://www.one.com/en/hosting/what-is-ssh?gclid=Cj0KCQiA1ZGcBhCoARIsAGQ0kkrTFwKizUcQFfXtTDS4WRrrzbZaVqN39hW5ROgMA8ilhZUerT15aWAaAnqlEALw_wcB)

To connect using ssh, we will use ssh-keys.

What are ssh-keys?(https://jumpcloud.com/blog/what-are-ssh-keys)

> **Exercise 1:** Generate a pair of ssh-keys. Send the public key to Ben for installation on the server along with your username. N.B. if working on Windows, try to use the Open SSH Client on the command line with the command: `ssh-keygen`. For both Windows and mac, the key should be saved to the default location which should be a hidden folder called `.ssh` inside your user's base directory.

> **Exercise 2:** ssh into the sequana server once Ben has installed your SSH keys. Ask Ben for sequana's IP. (optional) [set up](https://linuxize.com/post/using-the-ssh-config-file/) a 'config' file in '.ssh' so that you don't have to enter your IP and username everytime you want to 

## Part 2: Installing software in your user directory
We'll need some programmes if we're going to do some work.

One way to install programs on a Linux system is at the system level.

This allows all users to access the programs. Many of the programs you've already
used are installed at the systems level.

To see where a program is being executed from you can use the command `which`. E.g.:

```
$ which ls
/usr/bin/ls
```

While this might seem convenient there are many reasons not to install programs at the system level not limited to:
- You must have root access to be able to install the programs
- It is difficult to work with multiple versions of programs (i.e. v0.1.3 vs v0.1.4)

There are several ways to overcome the above issues. One is to download the source code for programs or pre-compiled programs and install these in a local directory. While this solves the problem of root access and the program will only be installed for your use, it still doesn't help us with managing multiple version of programs. It can also be difficult to download the source code of some programs and get it to compile correctly often due to dependency issues. However, sometimes downloading a program and installing it yourself is the only option available to you.

One of the first tasks we will need to undertake for our analysis is getting the sequencing data associated with the studies.

We will download the data from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home).

We will use their [File Downloader Command Line Tool](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html#using-ena-file-downloader-command-line-tool) that is hosted on GitHub [here](https://github.com/enasequence/ena-ftp-downloader/releases).

> **Exercise 3**: Install the ENA File Downloader in your home directory and make sure you can get it to execute.

We will work with this program later on to down load data. But first, let's meet another option for installing programs: [Conda](https://conda.io/projects/conda/en/latest/index.html) - a package manager that allows you to create multiple environments, often specific to a given project or analysis, that can have different packages and programs installed in them.

See [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) for managing environmens with conda. Conda can be installed in your home directory without root access (sudo access) and enables you to install programs locally.

> **Exercise 3**: Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) in your sequana home directory. If you're not sure which options to select, ask!

Great! We now have the ability to create environments and install programs locally using the 'conda' command. We will use conda later on.

> **Exercise 4**: Create a new environment called `testing_conda` and in it install: 
> - python3 v3.9
> - pandas (any version)
> - [kallisto](https://anaconda.org/bioconda/kallisto)
>
> Verify that you can run kallisto

## Part 3: Fetching data, directory permissions and working with symlinks
Now it's time to get some data.

> **Exercise 5**: Find where the data from each of the studies is deposited. Use the ENA File Downloader that you installed earlier to download data from ONE of the samples from the Böstrom et al 2017 paper.

Wow! Great work, you've now got real sequencing data in your home directory.

I only asked you to get 1 sample's worth of data because it would be redundant to have the same data downloaded multiple times on the server.

And guess what, I've already downloaded the data we need. It's in the following directory on sequana: XXX

One of the great things about operating on a Linux system is that its designed to allow concurrent access by multiple users.

> **Exercise 6**: Take a look at the permissions for the data directory. What does it tell you about access to the directory? What is a group in Linux? Which groups are you a member of?

Now you know where the data is hidden, you can either work directly with that data, or an easier way is to create a shortcut to that data in your own directory. For this purpose, Linux has symbolic links.

[What is a symbolic link (symlink)?](https://www.futurelearn.com/info/courses/linux-for-bioinformatics/0/steps/201767#:~:text=A%20symlink%20is%20a%20symbolic,directory%20in%20any%20file%20system.)

> **Exercise 7**: Create a directory structure in your home directory to hold the data. Create symlinks to populate the directories with symlinks to the sequencing files.

## Part 4: Some key sequence data formats - fasta, fastq, fastq.gz, sam and bam

There are a few key formats that you should be familiar with in the realms of computational biology.

[What is a fasta file?](https://en.wikipedia.org/wiki/FASTA_format)

[What is a fastq file?](https://en.wikipedia.org/wiki/FASTQ_format)

[What are sam and bam files?](https://www.zymoresearch.com/blogs/blog/what-are-sam-and-bam-files#:~:text=SAM%20files%20are%20a%20type,the%20examples%20for%20this%20section.)

# DAY 2: Böstrom et al 2017
## Part XXX:

# Day 3: MacParland et al 2018

## Requirements
You may want to download and install the [loupe browser](https://www.10xgenomics.com/products/loupe-browser/downloads). This is a browser designed by 10X to visualize the results of the Cell Ranger analysis that we will conduct below. Alternatively, I will load it on my machine and you can see what it looks like that way.

You'll be working with R to analyze the outputs of the scRNA-seq data analysis. If you'd like to run R on your local machine, you'll need to install R. You'll also likely want an IDE to work with. I would recommend either Rstudio or Visual Studio Code. You'll also need to install the following packages:

- dplyr
- Seurat
- patchwork
- ggplot2
- languageserver

Alternatively you can run R on the computation server (where it is already installed). But you'll likely want to connect to it using [Visual Studio Code](https://code.visualstudio.com/download). So install this on your local machine and install the R extension and the SSH extension.
## Part XXX: Running cellranger bamtofastq

TODO incorporate down sampling into this. If the students can down sample then they can work with their own sets of files and fully develop a Nextflow pipeline.

The processing of single cell RNA-seq data differs from that of bulk RNA-seq data with regards to the programs that are used to process it. However, many of the underlying principles (i.e. mapping and counting) are shared between bulk and single cell RNA-seq.

There are multiple ways to prepare single cell RNA-seq libraries. One popular way is to use the [Chromium controller from 10X Genomics](https://www.10xgenomics.com/instruments/chromium-family?utm_medium=search&utm_source=google&utm_campaign=sem-goog-2022-website-page-ra_g-chromium-brand-emea&useroffertype=website-page&userresearcharea=ra_g&userregion=emea&userrecipient=customer&usercampaignid=7011P000001mDXwQAM&gclid=Cj0KCQiA-JacBhC0ARIsAIxybyNyz7PZA_xjMRPViiRdmhkQQ2tuiga0crSjwBUrNtm5jJnddngIDbQaAlUiEALw_wcB&gclsrc=aw.ds).
This is what was used in MacParland et al 2018.

We will be working with a program called Cell Ranger which was developed by 10X Genomics.

[What is Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)?

The Cell Ranger pipeline generally starts from the raw output of the Illumina sequencing platform (raw base call format; BCL). Using the `cellranger mkfastq` command a set of fastq files are generated from the BCL. The set of fastq files is generated in a particular way where the various index sequences (i.e. barcode and unique molecular identifier (UMI)) are held in separate files to the RNA sequences. In Cell Ranger language, a barcode is unique to a cell (and a GEM well), and a UMI is unique to a transcript (pre-PCR).

Take a look at the data we downloaded from MacParland. You'll notice that it's in bam format, not fastq. This is normal for submission of 10X data to sequencing archives like NCBI. The bam format contains all the information of the fastq files, but with additional information on mapping.

As we touched on earlier, bam format is binary and not human readable. In order to read a bam file you will have to convert it to sam format.

> **Exercise XXX**: Use [samtools](http://www.htslib.org/) to look at the first 10 lines of the bam files in human readable format. How will you install it? Would you use conda or install it from source?

To run the Cell Ranger analysis of our data we want to run the `cellranger count` method. Take a look at the diagram in the below link to see what it is doing.

[What is `cellranger count` doing?](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview)

[Here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) is another good resource outlining the Single-Library Analysis that we'll be conducting.

As input, `cellranger count` takes fastq files. But we have bam files. Fortunately, Cell Ranger also contains a program to recapitulate the fastq files from bam files. It is helpfully named [bamtofastq](https://support.10xgenomics.com/docs/bamtofastq). We can run this now.

> **Excercise XXX**: Install cellranger. Find the bamtofastq executable and make sure you can get it to run.

We won't run it on the actual data because that would take a lot of time and resources that we don't have available to us right now. Instead I've recreated the output files here: `/home/humebc/VTK_22/macparland/results/bamtofastq`

## Part XXX: Running cellranger count
Now that we have the fastq files we can run the `cellranger count` program to generate the count tables.

> **Excercise XXX**: Run `cellranger count` on one set of the fastq files. You'll need to use a reference transcriptome. These come prebuilt for certain organisms (e.g. Human) and I have already downloaded it for you. It is here: `/home/humebc/VTK_22/reference/refdata-gex-GRCh38-2020-A`

An explanation of the output files can be found [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/gex-outputs)

> **Exercise XXX**: Examine the output of `cellranger count` and compare it to the above documentation.

## Part XXX: Getting files in and out of the server
Obviously, there will be many occasions where you want to either load files on to the server, or take files off of their server.

We will cover that now.

The most common way to do this is using [scp](https://linuxize.com/post/how-to-use-scp-command-to-securely-transfer-files/) (secure copy).

> **Exercise XXX**: Transfer the .html file that was generated to your local machines hard drive. Open it up.

> **Exercise XXX**: If you installed the Loupe browser, open up the Loupe file on your local machine. If you didn't, watch me!

## Part XXX: Running cellranger aggr

Great! Now we have one sample analyzed, but what about the other samples?

This is where `cellranger aggr` is used. [Here's](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate) the documentation on running aggr. And [here's](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/aggr-outputs) the documentation for the outputs.

> **Exercise XXX**: Do what's necessary to run `cellranger aggr`. Is it annoying having to run each sample individually on the command line? Imagine if you had 85 samples? Now how annoying is it? Do you think there are better approaches to running many samples, not just for Cell Ranger, but in general?

The documentation of the output structure for cellranger aggr is [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/aggr-outputs).


> **Exercise XXX**: Again, pull down the .html and cloupe files to visualize the results. Can you see the addition of the additional samples?

## Part XXX: Preparing R
We've done the heavy computing on the computational server. If you remember, the sequencing files that we started with were very large and the computation would have taken a very long time on your laptop - if it was possible at all.

However, the files that we've ended up with are not so large - about 100MB.

From here, the computations are less intense and the need for a large parallelization is reduced. We'll doing the next set of computations, and generating figures, in R.

[What is R?](https://www.r-project.org/about.html)

At this point we have a choice. We can either move the files that we generated and that are required for the next analysis onto our local machines and work with them there.

Alternatively, we can again 'connect' to the server and perform our analysis there.

If you choose to run on your local machine, then you'll need to install R and all of the dependencies that we'll be requiring. See the requirements sections at the beginning of this Day's section.

If you choose work on the server then you can either start up an R session on the command line (`R`) or you can use an IDE such as Visual Studio Code (recommended) to connect over SSH. This will make the work much easier for you.

> **Exercise XXX**: Get your chosen environment setup for the remainder of the analysis in R.

## Part XXX: scRNA-seq analysis with Seurat in R
Seurat is a package used to analyze scRNA-seq data. It has gained great popularity in recent years and is widely used by the academic community and industry alike.

Much of the work we'll be doing to recreate the results of the MacParland analysis are well documented by the creators of the Seurat package. For example, much of what we'll be doing is covered in their [pbmc3K tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).

> **Exercise XXX: Look through the methods of the MacParland paper to see where we're at with the analysis. Critically appraise how they've written up the methods. Is it easy to follow? Is there enough detail?

The first stage of analysis in R is to create a Seurat object from the features/barcode table that we created using the cellranger aggr command.

The table that we want to import is here: `/home/humebc/VTK_22/macparland/results/aggr/macparland/outs/count/filtered_feature_bc_matrix`

The commands for creating a Seurat object is as follow:
```
# First read in the data
mc.data <- Read10X(data.dir = "/home/humebc/VTK_22/macparland/results/aggr/macparland/outs/count/filtered_feature_bc_matrix", min.cells = 3)
```

Note that we are screening for features that are found in at least 3 cells, the same as they did in MacParland.

Until now we haven't really been able to get a feel for what the count table looks like. Now we can:

```
# Preview the table
mc.data[1:5, 1:30]

# Then make the Seurat object
mc <- CreateSeuratObject(counts = mc.data, project = "macparland")
```

You're welcome to call the Seurat variable whatever you like, but its probably easiest if you call it `mc` like I have.

A common part of processing the scRNA-seq data is to filter the data for 'high quality' cells.

There are several common paramaters by which cells can be filtered. E.g. see [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/).

One common parameter is the percentage of UMIs that are mitochondrial in origin.

We can calculate this metric for our data as follows:
```
mc[["percent.mt"]] <- PercentageFeatureSet(mc, pattern = "^MT-")
```

Other metrics have already been calculated for the dataset as part of creating the Seurat object:
```
head(pbmc@meta.data, 5)
summary(mc@meta.data)
```

It can be helpful to visualise these metrics:
```
# Create a violin plot
VlnPlot(mc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualise the correlation between metrics
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

Finally, we want to perform the filtering out cells from the data in the same way MacParland did. I.e.:
- filter out cell with < 1500 UMIs
- filter out cells with a high percentage of counts of mitochondrial origin

> **Exercise XXX**: Filter out those cells. Have a look at the pbmc3K tutorial for how to do this.

Now that we have completed the pre-processing of the data.

> **Exercise XXX**: Visualize the data again to make sure that the filtering has been applied. Has it?

> **Exercise XXX**: How would we produce the figure from the paper that shows library size plotted against mitochrondrial transcript percent?
## Part XXX: Normalization and dimensionality reduction

You'll hear the term normalization a lot in computational biology and data science.

It can have different meanings depending on the context.

Here, we're normalizing the feature expression measurements for each of the individual cells according to the total expression of the cell. E.g. think about a case where one cell has twice the number of UMIs sequenced than another another. If for a given feature the first cell has a count of 4 and in the second the cell has a count of 8, is that feature upregulated in the second cell?

The authors performed normalization using an R package called [scran](https://bioconductor.org/packages/release/bioc/html/scran.html). However, for simplicity's sake, we will continue using Seurat.

Normalize using Seurat:
```
mc <- NormalizeData(mc, normalization.method = "LogNormalize", scale.factor = 10000)
```

From here we are interested in clustering the data. That is assigning each of the cells to a given group.

To do this clustering we need to get an idea for how related each of the cells are to all of the other cells.

> **Exercise XXX**: Discuss, what are we basing this similarity on?

> **Exercise XXX**: Inspect the data, how many cells do we have? How many features?

Generally researchers choose to work on a subset of features for performing similarity analyses and clustering. This is to minimise computational load. All features could be used but do you think that all features provide the same ammount of information with regards to how similar cells are?

The most informative features will be those that vary the most between cells. In other words, they will be those features with the highest count variance across the cells.

We will identify these highly vairable features and use them for our down stream analyses:

```
mc <- FindVariableFeatures(mc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

The next step for clustering is to perform dimensionality reduction.

[What is dimensionality reduction?](https://en.wikipedia.org/wiki/Dimensionality_reduction)

We will do this by Principal Component Analysis (PCA)

[What is PCA?](https://en.wikipedia.org/wiki/Principal_component_analysis)

Before performing the PCA we need to standardize. This step gives equal weight to each of the features in downstream analyses, so that highly-expressed genes do not dominate.

Standardization shifts the expression of each gene, so that the mean expression across cells is 0 and scales the expression of each gene, so that the variance across cells is 1.

In Seurat this standardization is referred to as scaling.

Let's scale the data:
```
# Also referred to as standardization
all.genes <- rownames(mc)
mc <- ScaleData(mc, features = all.genes)
```

Then we can perform PCA
```
mc <- RunPCA(mc, features = VariableFeatures(object = mc))
```

And visualize the resultant Principal Components (PCs)
```
VizDimLoadings(mc, dims = 1:2, reduction = "pca")
ggsave("visdimloadings.png")

DimPlot(mc, reduction = "pca")
ggsave("dimplot.png")

# NB fast=F must be set in order to return a ggplot object.
DimHeatmap(mc, dims = 1, cells = 500, balanced = T, fast=F, nfeatures=100)
ggsave("dimheatmap.png")
```

Moving forward, we want to use a certain number of the resultant PCs. There wouldn't be much point in using all of the PCs. If we did that we wouldn't have reduced the dimensionaly of the dataset and the computational complexity would still be very high.

How many PCs do we select though?

One of the simplest ways to assess this is with an elbow plot:
```
ElbowPlot(mc, ndims=35)
```
> **Exercise XXX**: How many PCs should we move forwards with? How many did the authors choose?

## Part XXX: Clustering, non-linear dimensional reduction and identifying cluster biomarkers
Seurat implements a graph-based approach to clustering.

Cells are embeded in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.

The first part of this process is to find the Nearest Neibours of the cells:
```
mc <- FindNeighbors(mc, dims = 1:10)
```

Clusters are then called according to the Seurate algorithm from the graph:
```
mc <- FindClusters(mc, resolution = 0.5, pc.use=1:10)
```

Take a look at the cluster assignments:
```
head(Idents(mc), 5)
```

Finally we can run non-linear dimensional reduction in the form of UMAP or tSNE:
```
DimPlot(mc, reduction = "umap", pt.size=1)
ggsave("umap.png")

mc <- RunTSNE(mc, dims = 1:10)

DimPlot(mc, reduction = "tsne", pt.size=1)
ggsave("tsne.png")
```

In MacParland et al they identified the cells in the clusters by looking at cluster biomarkers. That is, genes that are unique markers of the cluster compared to all other clusters. They then used those identified genes with a manually curated set of genes that are kown to be indicative of certain cell types to assign cell types to the clusters.

We will not identify the cell type here, but we will identify the biomarkers (features) of the clusters.

Seurat includes power functionality to this end both to identify the markers and visually display the results for each of the clusters. We will produce an overview figure here as the last part of the practical:

```
mc.markers <- FindAllMarkers(mc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


mc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(mc, features = top10$gene) + NoLegend()
ggsave("cluster_diff_expression.png")
```

# Conclusions
Congratulations, you have completed the practical. I hope you have enjoyed it and found it informative. I will leave this practical online for your reference. Thanks for attending.













