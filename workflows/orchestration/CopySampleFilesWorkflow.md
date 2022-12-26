# Copy Sample Files Workflow
Workflow that copies sample files from one location to another.  See `FileUtils.CopySampleFiles`
task for more details.

### Input Parameters
* String sample_id - required - the id of the sample the files are for, passed to underlying task as a file match key
* String source_location - required - the source location to copy from
* String target_location - required - the target location to copy to
* Boolean flatten - optional - if true, don't replicate the relative directory structure in source location; defaults to false
* Boolean recursive - optional - if true, search recursively for files to copy; defaults to true
* String mgbpmbiofx_docker_image - required - the name/tag of the mgbpmbiofx/orchutils Docker image
* String gcp_project_id - required - the GCP project id to fetch location connection information from
* String workspace_name - required - the current Terra workspace name, used to fetch location connection information

### Output Parameters
* Array[String] source_files - a list of files copied from source
* Array[String] target_files - a list of files placed in target
* Array[File] local_files - if the target location is a local path (e.g. "./"), a list of the files
