process bam_to_junc {
    container = 'quay.io/karlvaba/leafcutter'

    input:
    file bam
    file bam_index

    output:
    path "${bam.baseName}.junc", emit: junc
    script:
    """
    regtools junctions extract -a 8 -m 50 -M 500000 $bam -o $bam.junc
    """
}

process intron_clustering_pr {
    container = 'quay.io/karlvaba/leafcutter'
    tag "${junc_files.baseName}"
    publishDir "${params.outdir}/leafcutter", mode: 'copy'

    input:
    file junc_files
    output:
    path "leafcutter_perind*.gz", emit: perind_counts
    script:
    """
    leafcutter_cluster.py -j $junc_files -m 50 -o leafcutter -l 500000 --checkchrom TRUE
    """
}


workflow intron_clustering {
    take:
    bam
    bam_index
    main:
    bam_to_junc(bam, bam_index)
    intron_clustering_pr(bam_to_junc.out.junc.map{it.toString()}.collectFile(name: 'junction_files.txt', newLine: true))
}