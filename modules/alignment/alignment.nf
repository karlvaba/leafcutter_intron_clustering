if (params.aligner == 'star') {
    include { star_align } from './star'
} else { //otherwise align with hisat2
    include { hisat2_align } from './hisat2'
}

workflow alignment {
    take:
        trimmed_reads
    main:
        if (params.aligner == 'star') {
            star_align(trimmed_reads)
        } else { //otherwise align with hisat2
            hisat2_align(trimmed_reads)
        }
    emit:
        bam = params.aligner == 'star' ? star_align.out.bam : hisat2_align.out.bam
        bam_index = params.aligner == 'star' ? star_align.out.bam_index : hisat2_align.out.bam_index
        
}