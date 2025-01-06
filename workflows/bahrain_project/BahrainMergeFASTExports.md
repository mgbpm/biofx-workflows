# Bahrain Merge Exports Workflow

The Bahrain Merge Exports Workflow will merge FAST exports for either the monogenic or screening pipeline and then run the FAST output parser.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | annotated_sample_data_name | Yes | Name for created annotated sample data in FAST | |
| Array[Files] | fast_export_files | Yes | Full paths to fast export files (result of FAST Parsing workflow) | |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| String | python_docker_image | No | The name of the Python Docker image | "python:3.10" |
| String | fast_parser_image | No | The name of the Docker image to run the FAST output parser task | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20240130" |
| File | portable_db_file | No | A SQLite database that contains additional annotations that are merged into the Parser output | "gs://lmm-reference-data/annotation/gil_lmm/gene_info.db" |
| String | fast_parser_sample_type | Yes | The sample type flag for the FAST output parser: S for single-sample Exome or M for multi-sample Exome or B for batch/Biobank or N for NVA-Lite | |
| Boolean | gatk_source | No | Whether or not to use GATK as a source for allele state annotations; true will use GATK as a source | false |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | merged_export_file | Always | Merged tab-delimited export of annotated sample data from FAST |
| File | fast_parsed_output | Always | Parsed FAST export |
