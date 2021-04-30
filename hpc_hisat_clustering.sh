../nextflow run ./main.nf\
 -profile tartu_hpc\
 --readPathsFile ./hpc_test_reads.tsv\
 --unstranded\
 --gtf /gpfs/hpc/projects/genomic_references/annotations/eQTLCatalogue/v0.1/Homo_sapiens.GRCh38.96.gtf\
 --hisat2_index /gpfs/hpc/projects/genomic_references/annotations/eQTLCatalogue/v0.1/hisat2_index_v96/Homo_sapiens.GRCh38.dna.primary_assembly\
 --aligner 'hisat'\
 -resume