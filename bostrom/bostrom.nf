#!/usr/bin/env nextflow

/*
All code associated with the analysis of the Böstrom dataset for the 2022 VTK
*/

samples_ch = Channel.fromFilePairs("/home/humebc/VTK_22/bostrom/raw_reads/*/*.fastq.gz", size: 1).map{[it[0], it[1][0]]}
grch38_index = file("/home/humebc/VTK_22/reference/kallisto_ensembl_ref/Homo_sapiens.GRCh38.cdna.all.short_name.kallisto.index")

process fastp{
    tag "${sample}"
    conda "fastp"
    publishDir "/home/humebc/VTK_22/bostrom/fastp/", pattern: "*.html"

    input:
    tuple val(sample), path(read_1) from samples_ch

    output:
    file "${sample}.fastp.html" into fastp_out_ch
    tuple val(sample), file("${sample}.clean.fq.gz") into kallisto_in_ch

    script:
    """
    fastp -q 20 -i $read_1 -o ${sample}.clean.fq.gz
    mv fastp.html ${sample}.fastp.html
    """
}

process kallisto{
    tag "${sample}"
    container "jennylsmith/kallistov45.0:latest"
    publishDir "/home/humebc/VTK_22/bostrom/kallisto/${sample}"
    cpus 10

    input:
    tuple val(sample), path(read_1_clean) from kallisto_in_ch
    file grch38_index

    output:
    file "*" into kallisto_out_ch

    script:
    """
    kallisto quant -i $grch38_index -o . -t ${task.cpus} --single -l 200 -s 30 $read_1_clean
    """
}