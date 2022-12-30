# Download Outputs Workflow
Wrapper around the `FileUtils.DownloadOutputsTask`.

### Input Parameters
* String outputs_json - the workflow outputs as a stringified JSON object where the keys are the output names
* String config_json - the output copy mapping configuration, see `biofx-workflow-configurations` for more details
* String? default_target_location - the default location to copy files to if there is not a more specific configuration
* String mgbpmbiofx_docker_image - required - the name/tag of the mgbpmbiofx/orchutils Docker image
* String gcp_project_id - required - the GCP project id to fetch location connection information from
* String workspace_name - required - the current Terra workspace name, used to fetch location connection information

### Output Parameters
* File outputs_manifest - the manifest of outputs and where they were copied to
