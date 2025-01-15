# Bahrain PGx/Risk Workflow

The Bahrain PGx/Risk Workflow is a part of the screening pipeline. It can be run at any time separately from the other parts of the pipeline. However, if run at the same time as the sample prep portion of the workflow, it can fail by reaching the egress limit on Terra.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[String] | sample_ids | Yes | List of sample IDs, e.g. SM-MPPD5; are included in the vcf and cram file names in the iCGD Clinical workspace | |
| Array[String] | collaborator_sample_ids  | Yes | List of sample IDs, e.g. D-981108334-BH-4022-S1-A | |
| Array[String] | data_bucket | Yes | List of vcf (and cram for screening pipeline) data locations for each sample from iCGD Clinical workspace | |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20240625" |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| File | ref_fasta | Yes | The genome reference FASTA file; required for running PGx and Risk with the screening pipeline | |
| File | ref_fasta_index | Yes | The genome reference FAST index file; required for running PGx and Risk with the screening pipeline | |
| File | ref_dict | No | The genome reference dict file; required for running PGx and Risk with the screening pipeline | |
| File | dbsnp_vcf | No | dbSNP VCF file that matches the genome reference; required for running PGx and Risk with the screening pipeline | |
| File | dbsnp_vcf_index | No | dbSNP VCF file index; required for running PGx and Risk with the screening pipeline | |
| String | pgx_test_code | No | Test code that defines which pharmacogenomics report to generate | "lmPGX-pnlC_L" |
| String | pgx_docker_image | No | The name of the Docker image to generate the pharmacogenomics report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:dev-v3" |
| File | pgx_workflow_fileset | No | Tar file containing the pharmacogenomics reference data to generate the report | "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_files-20220118.tar" |
| File | pgx_roi_bed | No | BED file that defines the genomic regions to include in the pharmacogenomics analysis | "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_genotyping-chr-20220118.bed" |
| String | risk_alleles_test_code | No | Test code that defines which risk alleles report to generate | "lmRISK-pnlB_L" |
| String | risk_alleles_docker_image | No | The name of the Docker image to generate the risk alleles report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20240129" |
| File | risk_alleles_workflow_fileset | No | Tar file containing the risk alleles reference data to generate the report | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar" |
| File | risk_alleles_roi_bed | No | BED file that defines the genomic regions to include in the risk alleles analysis | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | pgx_CPIC_report | When screening pipeline is run | CPIC pharmacogenomics report |
| Array[File] | pgx_FDA_report | When screening pipeline is run | FDA pharmacogenomics report |
| Array[File] | pgx_genotype_xlsx | When screening pipeline is run | Full list of pharmacogenomics genotypes in XLSX format |
| Array[File] | pgx_genotype_txt | When screening pipeline is run | Full list of pharmacogenomics genotypes in TSV format |
| Array[File] | risk_alleles_report | When screening pipeline is run | Risk alleles report |
| Array[File] | risk_alleles_genotype_xlsx | When screening pipeline is run | Full list of risk allele genotypes in XLSX format |
| Array[File] | risk_alleles_genotype_txt | When screening pipeline is run | Full list of risk allele genotypes in TSV format |
