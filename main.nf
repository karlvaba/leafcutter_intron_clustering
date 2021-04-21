#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//include { hisat2_align } from './modules/alignment/hisat2'
include { trim_galore } from './modules/trim_galore'
include { intron_clustering } from './modules/leafcutter/leafcutter'
include { star_align } from './modules/alignment/star'

workflow {
    if (params.skipAlignment) {
        bam = Channel.fromPath("${params.bamPath}/*.bam").ifEmpty{exit 1, "No bam files found"}
        bam_index = Channel.fromPath("${params.bamPath}/*.bam.bai").ifEmpty{exit 1, "No bam index files found"}
        sQTL_mapping(bam, bam_index)
    } else {
        trim_galore()
        //hisat2_align(trim_galore.out.trimmed_reads)
        star_align(trim_galore.out.trimmed_reads)
        //intron_clustering(hisat2_align.out.bam, hisat2_align.out.bam_index)
        intron_clustering(star_align.out.bam, star_align.out.bam_index)
    }
}