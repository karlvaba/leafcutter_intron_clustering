if( params.vcf ){
    Channel
        .fromPath("${params.vcf}*vcf.gz")
        .ifEmpty { exit 1, "VCF file(s) not found: ${params.vcf}" }
        .set { vcf }
} else {
    exit 1, "No VCF path specified!"
}

process modify_vcf_pr {
    container = 'quay.io/biocontainers/bcftools:1.12--h3f113a9_0'

    tag "${vcf.baseName} - modification"


    input:
    file vcf

    output:
    path "${vcf.baseName}.modified.vcf.gz", emit: vcf_modified

    shell:
    '''
    bcftools view -Oz -s `bcftools query -l !{vcf} | head -n 1` !{vcf} | gunzip | sed -e 's/[0-9]|[0-1]/0|1/g' | gzip -c > !{vcf.baseName}.modified.vcf.gz
    '''
}

workflow modify_vcf {
    main:
        modify_vcf_pr(vcf)
    emit:
        vcf_modified = modify_vcf_pr.out.vcf_modified
}