## Biobank RoR Data Sample Loading Workflow
The Biobank Return of Results pipeline is meant to filter Biobank datasets and create VCFs for individual samples within the datasets. This step of the pipeline is intended to be run following the Biobank RoR Individual Sample Preparation Workflow (BiobankRoRSamplePrepWorkflow). It annotate the individual sample VCFs from the sample prep workflow and load the annotations to FAST. 

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | VCF file for annotation and loading to FAST |
| String | sample_ID | Yes | ID of sample in input VCF |
| String | dataset | Yes | Name of Biobank dataset | |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20231129" |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| String | qceval_project_type | Yes | The type of rules to apply for the QC evaluation task, one of "WGS", "WGS_DRAGEN", "WES" or "NONE" | "WES" |
| String | qceval_docker_image | No | The name of the Docker image to run the QC evaluation task | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20231005" |
| String | queue_folder | No | Folder for Biobank RoR concurrency queue | "gs://biofx-task-queue/biobank-ror-1" | 
| Int | queue_limit | No | The maximum number of concurrent tasks to allow running in the queue | 200 |
| Int | queue_wait_limit_hrs | No | Number of hours that a sample can wait in the queue before failing | 15 |
| Int | fast_data_load_wait_interval_secs | No | The number of seconds in between checks when waiting for FAST data loads to complete | 120 |
| Int | fast_data_load_wait_max_intervals | No | The maximum number of checks to perform when waiting for FAST data loads to complete | 300 |
| Boolean | has_haploid_sites | No | If true, modify the VCF file headers prior to FAST load to work around lack of support Number=G fields and haploid sites | false |
| String | sample_data_load_config_name | No | The FAST load configuration name for the sample data VCF | "Sample_VCF_PPM_Eval" |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | qceval_vcf_gz | Always | Target VCF file annotated with QC Evaluation |
| String | fast_sample_data_name | Always | The name used when loading the sample to FAST |
| String | queue_entry_id | Always | The ID for the task running in the concurrency queue |