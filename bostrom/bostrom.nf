#!/usr/bin/env nextflow

/*
All code associated with the analysis of the MacParland dataset for the 2022 VTK
*/

samples_ch = Channel.fromFilePairs("/home/humebc/VTK_22/bostrom/raw_reads/*/*.fastq.gz", size: 1).map{[it[0], it[1][0]]}

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

// process kallisto{
//     tag "${sample}"
//     container "jennylsmith/kallistov45.0:latest"
//     publishDir "/home/humebc/projects/20220329_moya/spis_csm/kallisto_quant_results/${type}/${sample}"
//     cpus 10

//     input:
//     tuple val(type), val(sample), path(read_1_clean), path(read_2_clean) from kallisto_in_ch
//     file host_kallisto_index
//     file zooxs_kallisto_index

//     output:
//     file "*" into kallisto_out_ch

//     script:
//     if (type == "host")
//     """
//     kallisto quant -i $host_kallisto_index -o . -t ${task.cpus} $read_1_clean $read_2_clean
//     """
//     else if (type == "zooxs")
//     """
//     kallisto quant -i $zooxs_kallisto_index -o . -t ${task.cpus} $read_1_clean $read_2_clean
//     """
// }