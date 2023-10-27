## Loading Alamut and Gnomad Workflow Module
The Biobank Return of Results pipeline is meant to filter Biobank datasets and create VCFs for individual samples within the datasets. This step of the pipeline is intended to be run following the Biobank RoR Individual Sample Preparation Workflow (BiobankRoRSamplePrepWorkflow). It will batch annotate a prepared VCF from the sample prep workflow using Alamut and gnomAD, then load the annotations to FAST.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | batch_annotation_input | Yes | The file to use for Alamut and gnomAD batch annotation; collective VCF from Data Prep workflow | |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230921" |
| String | bcftools_docker_image | No | The name of the bcftools Docker image for VCF annotation | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17" |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| File | alamut_db | Yes | The database file for Alamut batch | |
| File | alamut_fields_tsv | No | The file that defines how the Alamut output is transformed back to a VCF | |
| String | alamut_db_name | No | The database name for the Alamut batch ini file | "alamut_db" |
| String | alamut_server | No | The server name for the Alamut batch ini file | "a-ht-na.interactive-biosoftware.com" |
| String | alamut_port | No | The server port for the Alamut batch ini file | "80" |
| String | alamut_user_secret_name | No | The GCP secret name that contains the user stanza for the Alamut batch ini file | "alamut-batch-ini-user" |
| Int | alamut_queue_limit | No | The maximum number of concurrent Alamut batch processes permitted | 4 |
| String | alamut_queue_folder | No | The shared storage location for Alamut concurrency management | "gs://biofx-task-queue/alamut" |
| String | alamut_docker_image | No | The name of the Alamut Docker image using to run Alamut Batch task | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/alamut:20230630" |
| Boolean | alamut_save_working_files | No | Whether or not to retain intermediate Alamut Batch task files | false |
| File | gnomad_coverage_file | Yes | The gnomad coverage data file | |
| File | gnomad_coverage_file_idx | Yes | The gnomad coverage data index file | |
| Array[String] | gnomad_headers | Yes | List of VCF headers to add when annotating VCF with gnomad coverage data.  For example, ##INFO=<ID=DP_gnomadG,Number=1,Type=Float,Description="Read depth of GnomAD Genome"> | |
| String | gnomad_column_list | Yes | The column list to pass to bcftools annotate for gnomad coverage annotation | "CHROM,POS,INFO/DP_gnomadG" |
| Boolean | has_haploid_sites | No | If true, modify the VCF file headers prior to FAST load to work around lack of support Number=G fields and haploid sites | false |
| String | gnomad_data_load_config_name | Yes | The FAST load configuration name for the gnomad coverage VCF | |
| String | alamut_data_load_config_name | No | The FAST load configuration name for the Alamut annotated VCF | "Alamut" |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | alamut_vcf_gz | Always | Target VCF file with Alamut annotations |
| File | gnomad_vcf_gz | Always | Target VCF annotated with gnomad coverage data |