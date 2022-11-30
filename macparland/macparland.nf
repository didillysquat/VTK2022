#!/usr/bin/env nextflow

/*
All code associated with the analysis of the MacParland dataset for the 2022 VTK
*/

// TODO add in the container for cell ranger and run from there rather than being
// reliant on the local version of cell ranger.


sample_bams_ch = Channel.fromFilePairs("/home/humebc/VTK_22/macparland/raw_reads/**/*.bam{,.bai}")
chromium_transcript_ref = file("/home/humebc/VTK_22/reference/refdata-gex-GRCh38-2020-A")

// The authors have provided us with bam files rather than the original fastq files.
// To get back to fastq files so that we can re run the analysis
// we will run CellRangers bamtofastq
// container "m0zack/bamtofastq:0.1"

//NB the continer "m0zack/bamtofastq:0.1" gave us errors about eof
// the container nfcore/cellranger doesn't contain the bamtofastq
// referenced in the path but it is in the image at /opt/cellranger-7.0.0/bamtofastq
process bamtofastq{
    tag "${sample}"
    
    publishDir "results/bamtofastq/${sample}/"
    cpus 20

    input:
    tuple val(sample), path(bams) from sample_bams_ch

    output:
    tuple val(sample), path("fastqs/*/*.fastq.gz") into fastqs_ch

    script:
    """
    /home/humebc/VTK_22/macparland/cellranger/cellranger-7.0.1/lib/bin/bamtofastq --nthreads=${task.cpus} ${bams[0]} ./fastqs
    """
}

// From here we end up with directories which contain the collection of fastq files
// that can then be used in the cellranger count function.
process count{
    tag "${sample}"
    cpus 35
    memory '150GB'
    publishDir "results/count/"

    input:
    path chromium_transcript_ref
    tuple val(sample), path(fastqs) from fastqs_ch

    output:
    path "${sample}" into count_out_ch
    
    script:
    """
    cellranger count --id $sample --transcriptome $chromium_transcript_ref --fastqs ./ --localcores ${task.cpus} --localmem 150
    """
}

// From here we need to run aggr to merge all of the individual cellranger count outputs into 
// a single barcode/feature table

// For this we first need to create an input csv table to provide to aggr
process aggr{
    tag "aggr"
    cpus 200
    publishDir "results/aggr"
    container "nfcore/cellranger:7.0.0"

    input:
    path count_dirs from count_out_ch.collect()

    output:
    path "data.csv" into aggr_out_ch
    
    script:
    """
    echo "sample_id,molecule_h5" > data.csv
    find -L ~+ -type f -name "molecule_info.h5" | awk -F '/' '{print \$(NF-2)","\$0}' >> data.csv
    cellranger aggr --id macparland --csv data.csv --localcores ${task.cpus} --localmem 800
    """
}
