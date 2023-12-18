## Biobank RoR FAST Parsing Workflow
The Biobank Return of Results pipeline is meant to filter Biobank datasets and create VCFs for individual samples within the datasets. This step of the pipeline is intended to be run following the Biobank RoR FAST Parsing Workflow (BiobankRoRFASTParsing). It will merge multiple exported FAST output files before running it through the FAST output parser.

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | annotated_sample_data_name | Yes | Name for created annotated sample data in FAST | |
| Array[Files] | fast_export_files | Yes | Full paths to fast export files (result of FAST Parsing workflow) | |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| String | python_docker_image | No | The name of the Python Docker image | "python:3.10" |
| String | fast_parser_image | No | The name of the Docker image to run the FAST output parser task | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20231206" |
| File | gil_transcript_exon_count | No | A tab delimited file of transcript id and exon count | "gs://lmm-reference-data/annotation/gil_lmm/transcript_exonNum.txt" |
| String | fast_parser_sample_type | The sample type flag for the FAST output parser: S for single-sample Exome or M for multi-sample Exome or B for batch/Biobank or N for NVA-Lite | B |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | merged_export_file | Always | Merged tab-delimited export of annotated sample data from FAST |
| File | fast_parsed_output | Always | Parsed FAST export |