## Biobank RoR FAST Parsing Workflow
The Biobank Return of Results pipeline is meant to filter Biobank datasets and create VCFs for individual samples within the datasets. This step of the pipeline is intended to be run following the Biobank RoR Data Sample Loading Workflow (BiobankRoRSampleLoadingWorkflow). It will create annotated sample data, export the data, and parse it.

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | gcp_project_id | Yes | The GCP project to fetch secrets from | |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | Yes | The name of the orchestration utils Docker image for FAST and file movement tasks | |
| File | batch_sample_ids | Yes | A list of sample IDs in the batch | |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| Array[String] | fast_annotated_sample_data_regions | No | The list of regions to include in the FAST annotated sample data; each element is a "name:applyMask" pair | |
| Array[String] | fast_annotated_sample_data_scripts | No | The list of custom scripts to run on the FAST annotated sample data after creation | |
| String | fast_annotated_sample_data_saved_filter_name | No | The saved filter to apply to the FAST annotated sample data | |
| String | fast_parser_image | Yes | The name of the Docker image to run the FAST output parser task | |
| File | gil_transcript_exon_count | Yes | A tab delimited file of transcript id and exon count |
| String | fast_parser_sample_type | Yes | The sample type flag for the FAST output parser: S for single-sample Exome or M for multi-sample Exome or B for batch/Biobank or N for NVA-Lite |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | fast_export_file | Always | Tab-delimited export of annotated sample data from FAST |
| File | fast_parsed_output | Always | Parsed FAST export |
| File | nva_report | Always | NVA report Excel document |