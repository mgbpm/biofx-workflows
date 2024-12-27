## Biobank RoR FAST Parsing Workflow
The Biobank Return of Results pipeline is meant to filter Biobank datasets and create VCFs for individual samples within the datasets. This step of the pipeline is intended to be run following the Biobank RoR Data Sample Loading Workflow (BiobankRoRSampleLoadingWorkflow). It will create annotated sample data, export the data, and parse it.

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | annotated_sample_data_name | Yes | Name for annotated sample data that will be created in FAST | |
| Array[String] | fast_data_sample_names | Yes | Names of samples in FAST for which to create annotated sample data
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20231129" |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| Array[String] | fast_annotated_sample_data_regions | No | The list of regions to include in the FAST annotated sample data; each element is a "name:applyMask" pair | |
| Array[String] | fast_annotated_sample_data_scripts | No | The list of custom scripts to run on the FAST annotated sample data after creation | |
| String | fast_annotated_sample_data_saved_filter_name | No | The saved filter to apply to the FAST annotated sample data | |
| String | fast_parser_image | No | The name of the Docker image to run the FAST output parser task | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20241226" |
| File | portable_db_file | No | A SQLite database that contains additional annotations that are merged into the Parser output | "gs://lmm-reference-data/annotation/gil_lmm/gene_info.db" |
| String | fast_parser_sample_type | No | The sample type flag for the FAST output parser: S for single-sample Exome or M for multi-sample Exome or B for batch/Biobank or N for NVA-Lite | B |
| Boolean | gatk_source | No | Whether or not to use GATK as a source for allele state annotations; true will use GATK as a source | false |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | fast_export_file | Always | Tab-delimited export of annotated sample data from FAST |
| File | fast_parsed_output | Always | Parsed FAST export |