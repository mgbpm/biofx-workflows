# Running a bamsurgeon workflow

## Background
<u>Usage</u>
bamsurgeon is a tool that will introduce genomic variants (such as SNV, SV, and indel mutations) into BAM, CRAM, and SAM files. The resulting mutated BAM files can be used to test variant callers.

<u>About the Workflow</u>
Bamsurgeon was originally developed by Adam Ewing. This workflow adapts the original Dockerfile and python scripts to work with a WDL on Terra. All scripts that are used in the tasks of this workflow can be found on the MGBPM-IT Azure DevOps Repo for bamsurgeon.

bamsurgeon has three main scripts, each corresponding to three categories of mutations that can be introduced to BAM files: addsnv.py, addindel.py, and addsv.py. Each of the different mutations has various optional inputs to run the workflow, most of which have defaults listed in the tables under the “Input Variables” section below.

Note that BWA makes assumptions that all genome reference files are conventionally named. The reference files found under the hg38 reference data in Terra suit this assumption. In addition, these reference files are not explicitly used in the WDL (besides the FASTA), but are necessary within the container for BWA to run.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[MutationBED] | mutation_bed | Yes | Array of target mutation information in JSON format | |
| File | basefile | Yes | BAM file to mutate | |
| File | basefile_idx | Yes | BAM index for BAM to mutate | |
| String | mutation_type | Yes | Type of mutation to introduce to input BAM file; either "snv", "sv", or "indel" | |
| String | output_basename | Yes | Base name for output files (should include sample name and mutation information) | |
| File | ref_fasta | Yes | Genome reference file | |
| File | ref_amb | Yes | Genome reference file | |
| File | ref_ann | Yes | Genome reference file | |
| File | ref_bwt | Yes | Genome reference file | |
| File | ref_fai | Yes | Genome reference file | |
| File | ref_pac | Yes | Genome reference file | |
| File | ref_sa | Yes | Genome reference file | |
| File | ref_dict | Yes | Genome reference file | |
| File | dbsnp_vcf | Yes | dbSNP VCF file that matches the genome reference | |
| File | dbsnp_vcf_index | Yes | dbSNP VCF file index | |
| Boolean | generate_reports | No | If `true`, the workflow will run PGx and Risk and generate the corresponding reports | false |
| String | pgx_test_code | No | Test code that defines which pharmacogenomics report to generate | "lmPGX-pnlD_L" |
| File | pgx_fileset | No | Tar file containing the pharmacogenomics reference data to generate the report | "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_20241004.tar" |
| File | pgx_roi | No | BED file that defines the genomic regions to include in the pharmacogenomics analysis | "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_genotyping.bed" |
| String | risk_test_code | No | Test code that defines which risk alleles report to generate | "lmRISK-pnlB_L" |
| File | risk_fileset | No | Tar file containing the risk alleles reference data to generate the report | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar" |
| File | risk_roi | No | BED file that defines the genomic regions to include in the risk alleles analysis | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed" |
| String | bamsurgeon_docker_image | No | Name of bamsurgeon Docker image | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bamsurgeon:20240229" |
| String | igvreport_docker_image | No | Name of IGV report Docker image | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511" |
| String | orchutils_docker_image | No | Name of orchutils Docker image | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20251222" |
| String | pgx_docker_image | No | The name of the Docker image to generate the pharmacogenomics report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20241007" |
| String | risk_alleles_docker_image | No | The name of the Docker image to generate the risk alleles report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20240129" |
| String | samtools_docker_image | No | Name of samtools Docker image | "biocontainers/samtools:v1.9-4-deb_cv1" |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | No | The name of the current workspace (for secret retrieval); require for running PGx/Risk | |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | mutated_bam | Always | Resulting sorted BAM file from bamsurgeon containing input SNV, indel, or SV mutation |
| File | mutated_bai | Always | Index file corresponding to the final sorted BAM file |
| File | igv_report | Always | IGV report detailing desired mutations in final BAM |
| File | pgx_summary_report | If `generate_reports` is true | PGx summarized Excel report with links to publications |
| File | pgx_details_report | If `generate_reports` is true | PGx detailed Excel report with Drug recommendations |
| File | pgx_genotype_xlsx | If `generate_reports` is true | Full list of pharmacogenomics genotypes in XLSX format |
| File | pgx_genotype_txt | If `generate_reports` is true | Full list of pharmacogenomics genotypes in TSV format |
| File | risk_report | If `generate_reports` is true | Risk alleles report |
| File | risk_genotype_xlsx | If `generate_reports` is true | Full list of risk allele genotypes in XLSX format |
| File | risk_genotype_txt | If `generate_reports` is true | Full list of risk allele genotypes in TSV format |
