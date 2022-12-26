# File Utility Tasks
The file utility tasks are wrappers around code implemented in the `biofx-orchestration-utils` repository.

## CopyFilesTask
Copies files from a source location to a target location.

### Input Parameters
* String source_location - required - the source location to copy from
* Array[String] file_types - optional - the list of file types to copy; if not specified, copy all
* Array[String] file_match_keys - optional - the list of strings that must be present in the path and name of each copied file
* String target_location - required - the target location to copy to
* Boolean flatten - optional - if true, don't replicate the relative directory structure in source location
* Boolean recursive - optional - if true, search recursively for files to copy
* String mgbpmbiofx_docker_image - required - the name/tag of the mgbpmbiofx/orchutils Docker image
* String gcp_project_id - required - the GCP project id to fetch location connection information from
* String workspace_name - required - the current Terra workspace name, used to fetch location connection information

### Output Parameters
* Array[String] source_files - a list of files copied from source
* Array[String] target_files - a list of files placed in target
* Array[File] local_files - if the target location is a local path (e.g. "./"), a list of the files; otherwise an empty array

## FetchFilesTask
Retrieves files from a location and places them in the task directory.

### Input Parameters
* String source_location - required - the source location to copy from
* Array[String] file_types - optional - the list of file types to copy; if not specified, copy all
* Array[String] file_match_keys - optional - the list of strings that must be present in the path and name of each copied file
* Boolean recursive - optional - if true, search recursively for files to copy
* String mgbpmbiofx_docker_image - required - the name/tag of the mgbpmbiofx/orchutils Docker image
* String gcp_project_id - required - the GCP project id to fetch location connection information from
* String workspace_name - required - the current Terra workspace name, used to fetch location connection information
* File empty_output_placeholder - optional - should never be specified; included only to support optional output parameters without runtime errors

### Output Parameters
* File? bam - the bam file, if present
* File? bai - the bam index file, if present
* File? cram - the cram file, if present
* File? crai - the cram index file, if present
* File? vcf - the vcf or gvcf file, if present
* File? vcf_index - the vcf or gvcf index file, if present
* Array[File] all_files - all retrieved files

