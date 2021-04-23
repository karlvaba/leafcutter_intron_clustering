#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { trim_galore } from './modules/trim_galore'
include { alignment } from './modules/alignment/alignment'
include { intron_clustering } from './modules/leafcutter/leafcutter'

workflow {
    trim_galore()
    alignment(trim_galore.out.trimmed_reads)
    intron_clustering(star_align.out.bam, star_align.out.bam_index)  
}