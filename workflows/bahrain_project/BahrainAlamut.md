# Bahrain Alamut Workflow

The Bahrain Alamut Workflow will take a merged or collective VCF from the the sample prep portion of the monogenic/screening pipeline and annotate with Alamut. It will then load the annotations to FAST.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | merged_vcf | Always | A merged VCF from the sample prep workflow | |
| String | batch_name | Yes | Prefix for all FAST sample data names, e.g. BGP-BH-4022 or BGP-Batch10 | |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20240625" |
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
| Boolean | has_haploid_sites | No | If true, modify the VCF file headers prior to FAST load to work around lack of support Number=G fields and haploid sites | false |
| String | alamut_data_load_config_name | No | The FAST load configuration name for the Alamut annotated VCF | "Alamut" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | alamut_vcf_gz | Always | Target VCF file with Alamut annotations |
| String | data_load_id | Always | Data load ID from loading Alamut annotations |
