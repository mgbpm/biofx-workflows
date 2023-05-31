# FAST Output Parser and NVA Report

# FASTOutputParserTask
The FASTOutputParserTask takes a tab-delimited export of FAST annotated sample data and produces
a parsed output with variants categorized that is used to populate an Excel workbook for review.

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | fast_output_file | Yes | The tab-delimited FAST annotated sample data export |
| String | sample_type | Yes | The sample type flag for the FAST output parser: S for single-sample Exome or M for multi-sample Exome or B for batch/Biobank or N for NVA-Lite |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| String | oms_query | Yes | Whether or not to query OMS for primers, either "Y" or "N" | "Y" |
| File | transcript_exonNum | Yes | A tab delimited file of transcript id and exon count |
| String | fast_parser_image | Yes | The name of the Docker image to run the FAST output parser task | |
| String | gcp_project_id | Yes | The GCP project to fetch secrets from | |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File? | parsed_report | Sometimes | The parsed FAST export with variants categorized or no file if FAST export contains no variants |
| File? | nva_report_xlsm | Sometimes | The output Excel workbook when the sample_type is not "N" and one or more variants are present |
| File? | nva_report_xlsx | Sometimes | The output Excel workbook when the sample_type is "N" or zero variants are present |
| File | nva_report | Always | The Excel workbook with the parsed report loaded |