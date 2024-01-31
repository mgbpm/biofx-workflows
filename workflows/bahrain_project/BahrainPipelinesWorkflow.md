# Bahrain Pipelines Workflow

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[String] | sample_ids | Yes | List of sample IDs, e.g. D-981108334-BH-4022-S1-A | |
| Array[String] | sample_data_locations | Yes | List of Wasabi data locations for each sample, e.g. s3://prod-bgp-1/BGP/Samples/BH-4022 | |
| String | batch_name | Yes | Prefix for all FAST sample data names, e.g. BGP_BH-4022 | |
| File | target_roi_bed | Yes | BED file containing gene regions; the VCFs in the dataset will be filtered to these regions and a new merged vcf will be created | |
| String | pipeline_to_run | Yes | Either "monogenic" to run the Bahrain Monogenic Pipeline or "screening" to run the Bahrain Screening Pipeline | |
| String | bcftools_docker_image | No | The name of the bcftools Docker image for VCF manipulation | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17" |
| String | ubuntu_docker_image | No | The name of the ubuntu Docker image | "ubuntu:latest" |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20231026" |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| File | ref_fasta | Yes | The genome reference FASTA file; required for running PGx and Risk with the screening pipeline | |
| File | ref_fasta_index | Yes | The genome reference FAST index file; required for running PGx and Risk with the screening pipeline | |
| File | ref_dict | No | The genome reference dict file; required for running PGx and Risk with the screening pipeline | |
| File | dbsnp_vcf | No | dbSNP VCF file that matches the genome reference; required for running PGx and Risk with the screening pipeline | |
| File | dbsnp_vcf_index | No | dbSNP VCF file index; required for running PGx and Risk with the screening pipeline | |
| String | pgx_test_code | No | Test code that defines which pharmacogenomics report to generate | "lmPGX-pnlC_L" |
| String | pgx_docker_image | No | The name of the Docker image to generate the pharmacogenomics report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20240129" |
| File | pgx_workflow_fileset | No | Tar file containing the pharmacogenomics reference data to generate the report | "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_files-20220118.tar" |
| File | pgx_roi_bed | No | BED file that defines the genomic regions to include in the pharmacogenomics analysis | "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_genotyping-chr-20220118.bed" |
| String | risk_alleles_test_code | No | Test code that defines which risk alleles report to generate | "lmRISK-pnlB_L" |
| String | risk_alleles_docker_image | No | The name of the Docker image to generate the risk alleles report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20240129" |
| File | risk_alleles_workflow_fileset | No | Tar file containing the risk alleles reference data to generate the report | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar" |
| File | risk_alleles_roi_bed | No | BED file that defines the genomic regions to include in the risk alleles analysis | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed" |
| File | alamut_db | Yes | The database file for Alamut batch | |
| File | alamut_fields_tsv | No | The file that defines how the Alamut output is transformed back to a VCF | |
| String | alamut_db_name | No | The database name for the Alamut batch ini file | "alamut_db" |
| String | alamut_server | No | The server name for the Alamut batch ini file | "a-ht-na.interactive-biosoftware.com" |
| String | alamut_port | No | The server port for the Alamut batch ini file | "80" |
| String | alamut_user_secret_name | No | The GCP secret name that contains the user stanza for the Alamut batch ini file | "alamut-batch-ini-user" |
| Int | alamut_queue_limit | No | The maximum number of concurrent Alamut batch processes permitted | 4 |
| String | alamut_queue_folder | No | The shared storage location for Alamut concurrency management | "gs://biofx-task-queue/alamut" |
| String | alamut_docker_image | No | The name of the Alamut Docker image using to run Alamut Batch task | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/alamut:20230630" |
| Boolean | alamut_save_working_files | No | Whether or not to retain intermediate Alamut Batch task files | false |
| String | alamut_anno_src_id | No | When removing already annotated variants prior to Alamut annotation, the annotation source id to query for | "228" |
| String | alamut_anno_min_age | No | When removing already annotated variants prior to Alamut annotation, the annotation minimum timestamp to query for (ISO8601 duration) | "P6M" |
| String | qceval_project_type | No | The type of rules to apply for the QC evaluation task, one of "WGS", "WGS_DRAGEN", "WES" or "NONE" | "WGS" |
| String | qceval_docker_image | No | The name of the Docker image to run the QC evaluation task | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20230511" |
| Boolean | has_haploid_sites | No | If true, modify the VCF file headers prior to FAST load to work around lack of support Number=G fields and haploid sites | false |
| String | sample_data_load_config_name | No | The FAST load configuration name for the sample data VCF | "Sample_VCF_PPM_Eval" |
| String | alamut_data_load_config_name | No | The FAST load configuration name for the Alamut annotated VCF | "Alamut" |
| Array[String] | fast_annotated_sample_data_regions | No | The list of regions to include in the FAST annotated sample data; each element is a "name:applyMask" pair | |
| Array[String] | fast_annotated_sample_data_scripts | No | The list of custom scripts to run on the FAST annotated sample data after creation | |
| String | fast_annotated_sample_data_saved_filter_name | No | The saved filter to apply to the FAST annotated sample data | |
| String | fast_parser_image | No | The name of the Docker image to run the FAST output parser task | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20240130" |
| File | gil_transcript_exon_count | Yes | A tab delimited file of transcript id and exon count |
| String | fast_parser_sample_type | No | The sample type flag for the FAST output parser: S for single-sample Exome or M for multi-sample Exome or B for batch/Biobank or N for NVA-Lite | "B" |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | prepped_sample_vcfs | Always | Sample VCFs that result from filtering to target ROI BED file and normalizing |
| File | merged_vcf_gz | Always | The resulting VCF from merging all prepped sample VCFs in a batch |
| File | collective_vcf_gz | When the screening pipeline is run | A VCF with a single "fake" sample but contains all the variants in the merged VCF |
| Array[File] | pgx_CPIC_report | When screening pipeline is run | CPIC pharmacogenomics report |
| Array[File] | pgx_FDA_report | When screening pipeline is run | FDA pharmacogenomics report |
| Array[File] | pgx_genotype_xlsx | When screening pipeline is run | Full list of pharmacogenomics genotypes in XLSX format |
| Array[File] | pgx_genotype_txt | When screening pipeline is run | Full list of pharmacogenomics genotypes in TSV format |
| Array[File] | risk_alleles_report | When screening pipeline is run | Risk alleles report |
| Array[File] | risk_alleles_genotype_xlsx | When screening pipeline is run | Full list of risk allele genotypes in XLSX format |
| Array[File] | risk_alleles_genotype_txt | When screening pipeline is run | Full list of risk allele genotypes in TSV format |
| File | alamut_vcf_gz | Always | Target VCF file with Alamut annotations |
| Array[File] | qceval_vcf_gz | Always | Target VCF file annotated with QC Evaluation |
| Array[File] | fast_export_file | Always | Tab-delimited export of annotated sample data from FAST |
| Array[File] | fast_parsed_output | Always | Parsed FAST export |
| Array[File] | nva_report | Always | NVA report Excel document |