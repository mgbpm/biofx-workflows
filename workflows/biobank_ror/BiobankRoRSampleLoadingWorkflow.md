## Biobank RoR Data Sample Loading Workflow
The Biobank Return of Results pipeline is meant to filter Biobank datasets and create VCFs for individual samples within the datasets. This step of the pipeline is intended to be run following the Biobank RoR Individual Sample Preparation Workflow (BiobankRoRSamplePrepWorkflow). It annotate the individual sample VCFs from the sample prep workflow and load the annotations to FAST. 

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | VCF file for annotation and loading to FAST |
| String | sample_ID | Yes | ID of sample in input VCF |
| String | dataset | Yes | Name of Biobank dataset | |
| String | batch | Yes | Name of batch or batch number | |
| String | gcp_project_id | Yes | The GCP project to fetch secrets from | |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | Yes | The name of the orchestration utils Docker image for FAST and file movement tasks | |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| String | qceval_project_type | No | The type of rules to apply for the QC evaluation task, one of "WGS", "WGS_DRAGEN", "WES" or "NONE" | "WES" |
| String | qceval_docker_image | Yes | The name of the Docker image to run the QC evaluation task | |
| Boolean | has_haploid_sites | No | If true, modify the VCF file headers prior to FAST load to work around lack of support Number=G fields and haploid sites | false |
| String | sample_data_load_config_name | Yes | The FAST load configuration name for the sample data VCF | |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | qceval_vcf_gz | Always | Target VCF file annotated with QC Evaluation |
| String | fast_sample_data_name | Always | The name used when loading the sample to FAST |