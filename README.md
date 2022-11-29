# Advanced Technologies for the Life Sciences: Bioinformatics Practical
The markdown document will cover all aspects of the bioinformatics
practical related to the 2022 VTK Advanced Technologies for the Life Sciences.

# Introduction
In this practical we will process 2 RNA sequencing (RNA-seq) datasets
to generate a selection of the analystical results presented in the two
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
The practical will take place in room M739 from 09.00-16.30 on the 6-8th of December. It is currently unclear whether we will have access to the room on the 9th of Decemeber as well.

Ben Hume, the current SequAna bioinformatician will
be running the practical and there to assisst you.

## Computing setup
You will need access to a computer to complete this practical.
While M739 does contain computers that you could use to login to the
SequAna computation server (where much of the work will occur)
this is a very inconvenient way to work as your home directory is
deleted after every session and it is not possible to install
applications that use Graphical User Interfaces (GUI).

As such it is strongly recommened that you bring a laptop with
you to complete the practical. There will be 1 spare laptop
available for use during the practical session in M739, but this
laptop cannot leave the room as it belongs to SequAna.

# Objectives

The main objective of this practcal is to give you an introduction to the tools that are used by computational biologist/bioinformaticians to generate meaningful results from sequencing data.

The objective of this course is not for you to become proficient or masterful of the techniques we will be covering (we have only 3 or 4 days!). Any proficiency gained, however, 
will likely be extremely valuable to you in your career as a research scientist.

To acheive this objective we will work with the sequencing data archived as part of the above mentioned stuies to recapitulate several of their key findings.

In doing so we will cover many broad informatic/bioinformatic techniques not limited to:

- Working on the command line interface (CLI)
- Using conda environments to install programs and packages
- Working with Docker images in Singularity
- Working with core bioinformatic tools to perform:
    - access of archived sequencing data
    - sequencing pre-processing and quality control
    - sequence analysis
- Workflow management with Nextflow
- R scripting to manipulate, analyse and visualise data

I will provide resources for all topics we cover and you are encouraged to look at these
resources if you wish to further your knowledge of the topic.

# Structure of the practical
The practical will be divided up by days (1-4).

Each day we will work towards our end goal of recapitulating the resutls of our chosen studies.

One of the most important skills in computation biology / informatics is the effective
sourcing of reference material. I.e. good googleing!

As such throughout the 4 days, while you will be given a sturcutre to follow,
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

This allows all users to access the progams. Many of the programs you've already
used are installed at the systems level.

To see where a program is being executed from you can use the command `which`. E.g.:

```
$ which ls
/usr/bin/ls
```

While this might seem convenient there are many reasons not to install programs at the system level not limited to:
- You must have root access to be able to install the programs
- It is difficult to work with multiple versions of programs (i.e. v0.1.3 vs v0.1.4)

There are several ways to overcome the above issues. One is to download the source code for programs or pre-compiled programs and install these in a local directory. While this solves the problem of root access and the program will only be installed for your use, it still doesn't help us with managing multiple version of programs. It can also be difficult to download the source code of some programs and get it to complile correctly often due to dependency issues. However, sometimes downloading a program and installing it yourself is the only option available to you.

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

TODO cover the formats

Can use conda to intall samtools to view the bam file to look at the structures.

# DAY 2: Böstrom et al 2017
## Part XXX:

# Day 3: MacParland et al 2018

TODO incorporate down sampling into this. If the students can down sample then they can work with their own sets of files and fully develop a Nextflow pipeline.

The processing of single cell RNA-seq data differs from that of bulk RNA-seq data with regards to the programs that are used to process it. However, many of the underlying principles (i.e. mapping and counting) are shared between bulk and single cell RNA-seq.

There are multiple ways to prepare single cell RNA-seq libraries. One popular way is to use the [Chromium controller from 10X Genomics](https://www.10xgenomics.com/instruments/chromium-family?utm_medium=search&utm_source=google&utm_campaign=sem-goog-2022-website-page-ra_g-chromium-brand-emea&useroffertype=website-page&userresearcharea=ra_g&userregion=emea&userrecipient=customer&usercampaignid=7011P000001mDXwQAM&gclid=Cj0KCQiA-JacBhC0ARIsAIxybyNyz7PZA_xjMRPViiRdmhkQQ2tuiga0crSjwBUrNtm5jJnddngIDbQaAlUiEALw_wcB&gclsrc=aw.ds).
This is what was used in MacParland et al 2018.

We will be working with a program called Cell Ranger which was developed by 10X Genomics.

[What is Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)?

The Cell Ranger pipeline generally starts from the raw output of the Illumina sequencing platform (raw base call format; BCL). Using the `cellranger mkfastq` command a set of fastq files are generated from the BCL. The set of fastq files is generated in a particular way where the various index sequences (i.e. barcode and unique molecular identifier (UMI)) are held in separate files to the RNA sequences. In Cell Ranger language, a barcode is unique to a cell (and a GEM well), and a UMI is unique to a transcript (pre-PCR).

Take a look at the data we downloaded from MacParland. You'll notice that it's in bam format, not fastq. This is normal for submission of 10X data. The bam format contains all the information of the fastq files, but with additional information on mapping.

As we touched on earlier, bam format is binary and not human readable. In order to read a bam file you will have to convert it to sam format.

> **Exercise XXX**: Use [samtools](http://www.htslib.org/) to look at the first 10 lines of the bam files in human readable format. How will you install it? Would you use conda or install it from source?

To run the Cell Ranger analysis of our data we want to run the `cellranger count` method. Take a look at the diagram in the below link to see what it is doing.

[What is `cellranger count` doing?](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview)

[Here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) is another good resource outlining the Single-Library Analysis that we'll be conducting.

As input, `cellranger count` takes fastq files. But we have bam files. Fortunately, Cell Ranger also contains a program to recapitulate the fastq files from bam files. It is helpfully named [bamtofastq](https://support.10xgenomics.com/docs/bamtofastq). We can run this now.

> **Excercise XXX**: Install cellranger. Find the bamtofastq executable and make sure you can get it to run.

We won't run it on the actual data because that would take a lot of time and resources that we don't have available to us right now. Instead I've recreated the output files here: XXX

Now that we have the fastq files we can run the `cellranger count` program to generate the count tables.

> **Excercise XXX**: Run `cellranger  count` on one set of the fastq files.





