# T1K HLA Typing Workflow

This workflow will take either a sample BAM or CRAM and run T1K HLA genotyping. Further documentation for T1K can be found on [Github](https://github.com/mourisl/T1K).

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_sample | Yes | Either a BAM or CRAM of sample data desired for HLA typing | |
| File | input_sample_idx | Yes | Respective CRAI or BAI file for the input_sample file | |
| File | ref_fasta | No | Reference GRCh38 FASTA file | "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta" |
| File | ref_fai | No | Respective index file for the reference GRCh38 FASTA file | "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai" |
| String | t1k_docker_image | No | Docker image with T1K software | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/t1k:0.0.1" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | allele_tsv | Always | Representative alleles with all fields and respective quality scores; To be used in conjunction with the allele_vcf output |
| File | allele_vcf | Always | Novel SNPs |
| File | genotype_tsv | Always | Additional metrics for alleles found in allele_tsv |
