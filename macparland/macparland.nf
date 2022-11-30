#!/usr/bin/env nextflow

/*
All code associated with the analysis of the MacParland dataset for the 2022 VTK
*/

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
    tuple val(sample), path("${sample}") into count_out_ch
    
    script:
    """
    cellranger count --id $sample --transcriptome $chromium_transcript_ref --fastqs ./ --localcores ${task.cpus} --localmem 150
    """
}

