if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .set { gtf }
} else {
    exit 1, "No GTF annotation specified!"
}

if ( params.fasta ){
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
           .set { fasta }
}
else {
    exit 1, "No reference genome specified!"
}

if ( params.star_index ){
    star_index = Channel
        .fromPath("${params.star_index}")
        .ifEmpty { exit 1, "Star index not found: ${params.star_index}" }
}

else {
    exit 1, "No star index specified!"
}



process makeSTARindex {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$fasta"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file fasta
    file gtf

    output:
    path "star", emit: star_index

    script:
    def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile $gtf \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        $avail_mem
    """
}

skipped_poor_alignment = []
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() <= '5'.toFloat() ){
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        skipped_poor_alignment << logname
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}

process star {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$prefix"
    publishDir "${params.outdir}/STAR", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    tuple val(samplename), file(reads)
    file star_index
    file gtf
    file vcf

    output:
 
    tuple file("*Log.final.out"), file ('*.bam'), emit: star_aligned
    path "*.out", emit: alignment_logs
    path "*SJ.out.tab"
    path "*Log.out", emit: star_log
    path "${prefix}Aligned.sortedByCoord.out.bam.bai", emit: bam_index

    script:
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    def star_mem = task.memory ?: params.star_memory ?: false
    def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''
    seqCenter = params.seqCenter ? "--outSAMattrRGline ID:$prefix 'CN:$params.seqCenter'" : ''
    if (params.wasp) {
        """
        STAR --genomeDir $star_index \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads  \\
            --runThreadN ${task.cpus} \\
            --twopassMode Basic \\
            --waspOutputMode SAMtag \\
            --outSAMattributes NH HI AS nM NM vA vG vW \\
            --varVCFfile <(zcat ${vcf}) \\
            --outWigType bedGraph \\
            --outSAMtype BAM SortedByCoordinate $avail_mem \\
            --readFilesCommand zcat \\
            --runDirPerm All_RWX \\
                --outFileNamePrefix $prefix $seqCenter
            
        samtools index ${prefix}Aligned.sortedByCoord.out.bam
        """
    } else {
        """
        STAR --genomeDir $star_index \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads  \\
            --runThreadN ${task.cpus} \\
            --twopassMode Basic \\
            --outWigType bedGraph \\
            --outSAMtype BAM SortedByCoordinate $avail_mem \\
            --readFilesCommand zcat \\
            --runDirPerm All_RWX \\
                --outFileNamePrefix $prefix $seqCenter
            
        samtools index ${prefix}Aligned.sortedByCoord.out.bam
        """
    }
    
    }

include { modify_vcf } from '../utils/vcf_mod'

workflow star_align {
    take:
        trimmed_reads
    main:
        modify_vcf()
        star(trimmed_reads, star_index, gtf.collect(), modify_vcf.out.vcf_modified)
    emit:
        bam = star.out.star_aligned.filter { logs, bams -> check_log(logs) }.flatMap {  logs, bams -> bams }
        bam_index = star.out.bam_index
}