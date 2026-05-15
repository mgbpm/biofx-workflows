# GDH Ingest and Filter Tasks

# GDHIngestAndFilterTask
Ingests a VCF file into the GDH (Genomic Data Hub) platform and runs a filter against the ingested data, returning matching variants.

## Input Parameters
* String subject_id - required - the subject/accession identifier
* String sample_id - required - the sample identifier
* String gdh_institution - optional - the institution name, defaults to "MGBPM"
* String gdh_project - optional - the project name, defaults to "Clinical"
* File vcf_file - required - the VCF file to ingest
* String vcf_file_stage_name - optional - the stage name for the VCF file, defaults to "biofx_pipelines"
* String vcf_file_stage_gspath - optional - the GCS path for the VCF file stage, defaults to "gs://gdh-external-stage/biofx_pipelines_nonprod"
* String reference_build - optional - the reference genome build, defaults to "GRCh38"
* Array[String] vcf_transform_functions - optional - list of VCF transform functions to apply, defaults to ["sample_base.lmm_calculate_variant_call_attributes"]
* String filter_name_or_code - required - the name or code of the filter to apply
* String? pipeline_run_id - optional - an identifier for the pipeline run; if not provided, defaults to subject_id + sample_id
* Int timeout_minutes - optional - the maximum time in minutes to wait for the ingest and filter to complete, defaults to 90
* String gcp_project_id - required - the GCP project to fetch secrets from
* String workspace_name - required - the name of the current workspace (for secret retrieval)
* String docker_image - required - the Docker image name:tag
* Int disk_size - optional - the amount of disk space to request, defaults to the VCF file size plus 10
* Boolean verbose - optional - enable verbose output, defaults to false

## Output Parameters
* String ingest_and_filter_execution_id - the execution identifier for the ingest and filter run
* File matching_variants - a JSON file containing the variants matching the specified filter

# GDHOutputParserTask
Parses GDH filter output into human-readable report formats (text and Excel).

## Input Parameters
* File gdh_output_file - required - the GDH output file to parse (JSON or text format)
* String reference_build - required - the reference genome build
* String oms_query - optional - whether to query OMS for additional information, defaults to "Y"
* File portable_db_file - required - the portable database file containing gene information
* String report_basename - optional - the basename for the output report files, defaults to the input file name with extensions removed
* String filter_name_or_code - required - the name or code of the filter that was applied
* String gcp_project_id - required - the GCP project to fetch secrets from
* String workspace_name - required - the name of the current workspace (for secret retrieval)
* String gdh_parser_image - required - the Docker image name:tag for the GDH parser

## Output Parameters
* File? parsed_report - the parsed text report
* File? nva_report_xlsx - the Excel report in .xlsx format
* File? nva_report_xlsm - the Excel report in .xlsm format (with macros)
* File nva_report - the preferred Excel report (xlsm if available, otherwise xlsx)
