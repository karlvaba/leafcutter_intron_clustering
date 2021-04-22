../nextflow run ./main.nf\
 -profile tartu_hpc\
 --readPathsFile ./hpc_test_reads.tsv\
 --reverse_stranded\
 --gtf /gpfs/hpc/projects/genomic_references/annotations/eQTLCatalogue/v0.1/Homo_sapiens.GRCh38.96.gtf\
 --fasta /gpfs/hpc/projects/genomic_references/annotations/eQTLCatalogue/v0.1/Homo_sapiens.GRCh38.dna.primary_assembly.fa\
 --wasp\
 -resume