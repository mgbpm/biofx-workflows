# Download Outputs Workflow
Workflow that distributes workflow output files to various locations as defined by the provided output mapping
configuration.  Intended to be used to download outputs for entity set workflows.  Only the first element of the
default target location list is used; the others are ignored.  The parameter accepts an array to be compatible
with a `this.[setmembers].default_target_location` expression in Terra.

See also `FileUtils.DownloadOutputsTask`.

### Input Parameters
* String outputs_json - the workflow outputs as a stringified JSON object where the keys are the output names
* Array[String] config_json_list - a list of output copy mapping configurations, see `biofx-workflow-configurations` for more details
* Array[String] default_target_location_list - list of default locations to copy files to if there is not a more specific configuration; defaults to empty array
* String orchutils_docker_image - required - the name/tag of the mgbpmbiofx/orchutils Docker image
* String gcp_project_id - required - the GCP project id to fetch location connection information from
* String workspace_name - required - the current Terra workspace name, used to fetch location connection information

### Output Parameters
* File outputs_manifest - the manifest of outputs and where they were copied to

