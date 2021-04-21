if(params.readPathsFile){
    if(params.singleEnd){
        Channel.fromPath(params.readPathsFile)
        .ifEmpty { error "Cannot find any readPathsFile file in: ${params.readPathsFile}" }
        .splitCsv(header: false, sep: '\t', strip: true)
        .map{row -> [ row[0], [ file(row[1]) ] ]}
        .set { raw_reads }
    } else {
        Channel.fromPath(params.readPathsFile)
        .ifEmpty { error "Cannot find any readPathsFile file in: ${params.readPathsFile}" }
        .splitCsv(header: false, sep: '\t', strip: true)
        .map{row -> [ row[0], [ file(row[1]) , file(row[2]) ] ]}
        .set { raw_reads }
    }
} 

process trim_galore_pr {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else null
        }

    input:
    tuple val(name), file(reads)

    output:
    tuple val(name), file("*fq.gz"), emit: trimmed_reads
    path "*trimming_report.txt", emit: trimgalore_results
    path "*_fastqc.{zip,html}", emit: trimgalore_fastqc_reports


    script:
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    if (params.singleEnd) {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}

workflow trim_galore {
    main:
        trim_galore_pr(raw_reads)
    emit:
        trimmed_reads = trim_galore_pr.out.trimmed_reads
}