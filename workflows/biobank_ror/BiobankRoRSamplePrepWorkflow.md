## Biobank RoR Individual Sample Preparation Workflow
The Biobank Return of Results pipeline is meant to filter Biobank datasets and create VCFs for individual samples within the datasets. This step of the pipeline will prep individual sample VCFs and create a TSV file to upload to the Biobank RoR Terra workspace for the next steps in the pipeline.

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[String] | filenames | Yes | List of input file names | |
| Array[String] | file_locations | Yes | List of Wasabi or GCP locations for input files | |
| String | fetch_file_type | Yes | Type of input file that is being fetched; either "vcf" or "bcf" | |
| File | sample_ids_list | No | A tsv file containing sample IDs from the dataset; if not provided, one will be created from the fetched files | |
| String | dataset | Yes | Name of Biobank dataset | |
| String | dataset_structure | Yes | Whether the dataset is made of joint VCFs or individual sample VCFS; either "individual" or "joint" | |
| File | target_roi_bed | No | BED file containing gene regions; the VCFs in the dataset will be filtered to these regions | |
| Boolean | is_sharded | Yes | Whether or not the input data is sharded | |
| String | python_docker_image | No | The name of the python Docker image for creating a collective VCF | "python:3.10" |
| String | bcftools_docker_image | No | The name of the bcftools Docker image for VCF annotation | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17" |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20231129" |
| String | ubuntu_docker_image | No | The name of the ubuntu Docker image | "ubuntu:latest" |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | filtered_vcfs | When filtering input VCFs is applicable | Result of filtering dataset VCFs |
| File | concat_vcf | When concatenating dataset VCFs is applicable; when the dataset is sharded | Result of concatenating dataset VCFs |
| Array[File] | individual_vcfs | When the dataset contains joint VCFs | Result of splitting the dataset VCFs into per sample VCFs |
| File | batch_annotation_input_file | Always | Collective VCF file to use for batch annotation with Alamut and gnomAD; used as input for BiobankRoRAlamutWorkflow |
| File | dataset_sample_table | Always | TSV file containing file paths for individual sample VCFs; upload to Terra to continue the pipeline |