# File Utility Tasks
The file utility tasks are wrappers around code implemented in the `biofx-orchestration-utils` repository.
See code documentation for more detailed behavior explanations.

## FileMatcher Struct
Defines a single matching rule for file selection.  Defines the following properties:
* file_type - the file extension
* path_matches - a list of strings or regular expressions that must be present in the relative file path
* required - whether a matching file must be found

## CopyFilesTask
Copies files from a source location to a target location.

### Input Parameters
* String source_location - required - the source location to copy from
* Array[String] file_types - optional - the list of file types to copy; if not specified, copy all
* Array[String] file_match_keys - optional - the list of strings that must be present in the path and name of each copied file
* Array[FileMatcher] file_matchers - optional - the list of file matching rules; mutually exclusive with file_types and file_match_keys
* String target_location - optional - the target location to copy to; defaults to the local working directory
* Boolean flatten - optional - if true, don't replicate the relative directory structure in source location
* Boolean recursive - optional - if true, search recursively for files to copy
* Boolean verbose - optional - if true, generate verbose log output
* String docker_image - required - the name/tag of the orchutils Docker image
* Int disk_size - optional - size of disk to allocation in GB, defaults to 75
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
* Array[FileMatcher] file_matchers - optional - the list of file matching rules; mutually exclusive with file_types and file_match_keys
* Boolean recursive - optional - if true, search recursively for files to copy
* Boolean verbose - optional - if true, generate verbose log output
* String docker_image - required - the name/tag of the orchutils Docker image
* Int disk_size - optional - size of disk to allocation in GB, defaults to 75
* String gcp_project_id - required - the GCP project id to fetch location connection information from
* String workspace_name - required - the current Terra workspace name, used to fetch location connection information
* File empty_output_placeholder - optional - should never be specified; included only to support optional output parameters without runtime errors

### Output Parameters
* File? bam - the bam file, if present; if multiple match, the first when sorted alphanumerically
* File? bai - the bam index file, if present; if multiple match, the first when sorted alphanumerically
* File? cram - the cram file, if present; if multiple match, the first when sorted alphanumerically
* File? crai - the cram index file, if present; if multiple match, the first when sorted alphanumerically
* File? vcf - the vcf or gvcf file, if present; if multiple match, the first when sorted alphanumerically
* File? vcf_index - the vcf or gvcf index file, if present; if multiple match, the first when sorted alphanumerically
* File? bcf - the bcf or bcfgz file, if present; if multiple match, the first when sorted alphanumerically
* Array[File] all_files - all retrieved files

## DownloadOutputsTask
Downloads outputs from a workflow using the provided mapping

### Input Parameters
* String outputs_json - the workflow outputs as a stringified JSON object where the keys are the output names
* String config_json - the output copy mapping configuration, see `biofx-workflow-configurations` for more details
* String? default_target_location - the default location to copy files to if there is not a more specific configuration
* Boolean verbose - optional - if true, generate verbose log output
* String docker_image - required - the name/tag of the orchutils Docker image
* String gcp_project_id - required - the GCP project id to fetch location connection information from
* String workspace_namespace - optional - the current Terra workspace namespace, used to fetch submission metadata
* String workspace_name - required - the current Terra workspace name, used to fetch location connection information and submission metadata
* String submission_id - optional - the submission id to fetch metadata and logs for

### Output Parameters
* File outputs_manifest - the manifest of outputs and where they were copied to

