# FAST Utility Tasks

# FASTDataLoadTask
Starts a sample or annotation data load using the specified parameters.

## Input Parameters
* File vcf_file - required - the VCF file to load
* String? default_annotation_src - optional - the name of the default annotation source
* String? reference_build - optional - the reference genome build
* String data_load_target - required - either SAMPLE_DATA or ANNOTATION_DATA
* String? sample_data_exists_strategy - optional - what to do when sample data already exists, defaults to FAIL
* String? merge_strategy - optional - rules to apply to combining annotations from for the same concept and source, defaults to MERGE
* String? sample_data_coll - optional - the sample data collection name
* String? sample_data_name - optional - the sample data name
* String? lab_batch_name - optional - the lab batch name
* String data_load_config_name - optional - the name of the configuration to use for data loading
* String? custom_script - optional - a custom script to execute following the data load
* String? annotation_record_ts - optional - fixed timestamp to set for annotation records, in ISO8601 format
* String? email_to - optional - email to send completion notice to
* String gcp_project_id - required - the GCP project to fetch secrets from
* String workspace_name - required - the name of the current workspace (for secret retrieval)
* String docker_image - required - the Docker image name:tag
* Int disk_size - optional - the amount of disk space to request, defaults to 20

## Output Parameters
* String data_load_id - the identifier for the data load created

# FASTWaitForDataLoadsTask
Waits up to a certain amount of time for a data load to complete.

## Input Parameters
* Array[String] data_load_ids - required - list of data load ids to wait for
* Int check_interval_secs - optional - how frequently (in seconds) to check data load status, defaults to 60
* Int max_checks - optional - how many times to check status before giving up, defaults to 24 * 60
* String gcp_project_id - required - the GCP project to fetch secrets from
* String workspace_name - required - the name of the current workspace (for secret retrieval)
* String docker_image - required - the Docker image name:tag

## Output Parameters
* FASTWaitResult wait_result - counts of the number of data loads: total, completed, succeeded, failed, canceled

# FASTCreateAnnotatedSampleDataTask
Initiates the creation of an annotated sample data

## Input Parameters
* String annotated_sample_data_name - required - name of the annotated sample data
* Array[String] sample_data_names_and_labels - required - list of sample data names to include in the annotated sample data;
  for each sample data, a label may be added by appending it to the sample data name, e.g. "Sample:Label"
* Array[String]? region_names_and_masks - optional - a list of target region names to include in the annotated sample data;
  for each region, whether to apply the region as a mask may be specified by appending it to the region name, e.g. "Region:true";
  if a mask value is not specified, it defaults to false
* Array[String]? scripts - optional - a list of custom scripts to execute on annotated sample data creation
* String? saved_filter_name - optional - the saved filter to apply to the annotated sample data
* String gcp_project_id - required - the GCP project to fetch secrets from
* String workspace_name - required - the name of the current workspace (for secret retrieval)
* String docker_image - required - the Docker image name:tag
* Int disk_size - optional - disk size for the container in GB, defaults to 20

## Output Parameters
* String annotated_sample_data_name_output - a copy of the input parameter, available to force task execution dependency

# FASTWaitForAnnotatedSampleDataTask
Waits for completion of an annotated sample data

## Input Parameters
* String annotated_sample_data_name - required - name of the annotated sample data
* Int check_interval_secs - optional - how frequently (in seconds) to check data load status, defaults to 60
* Int max_checks - optional - how many times to check status before giving up, defaults to 24 * 60
* String gcp_project_id - required - the GCP project to fetch secrets from
* String workspace_name - required - the name of the current workspace (for secret retrieval)
* String docker_image - required - the Docker image name:tag

## Output Parameters
* FASTWaitResult wait_result - counts of the number of data loads: total, completed, succeeded, failed, canceled

# FASTRemoveAlreadyAnnotatedFromVCFTask
Removes rows from a VCF file where the annotation data in FAST already includes at least one annotation matching
the specified criteria.  At least one annotation criterion must be specified.

## Input Parameters
* File input_vcf - required - the input VCF file, may be compressed
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
  VCF name with file extensions removed and ".alreadyann" added
* String reference_build - required - the reference build to use when querying FAST
* String annotation_source_id - optional - the annotation record source id
* String annotation_field_id - optional - the annotation record field id
* String annotation_min_age - optional - the annotation record minimum timestamp, ISO8601 timestamp or duration
* String gcp_project_id - required - the GCP project to fetch secrets from
* String workspace_name - required - the name of the current workspace (for secret retrieval)
* String docker_image - required - the Docker image name:tag
* Int batch_size - optional - the number of rows to process in each pass (FAST query), defaults to 20
 
## Output Parameters
* File output_vcf_gz - the output VCF file

## Docker Image Requirements
* Python3
* biofx-orchestration-utils code repository
* biofx-pyfast code repository