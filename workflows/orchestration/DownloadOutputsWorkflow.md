# Download Outputs Workflow
Workflow that distributes workflow output files to various locations as defined by the provided output mapping
configuration.  Intended to be used to download outputs for single entity workflows.

See also `FileUtils.DownloadOutputsTask`.

### Input Parameters
* String outputs_json - the workflow outputs as a stringified JSON object where the keys are the output names
* String config_json - the output copy mapping configuration, see `biofx-workflow-configurations` for more details
* String default_target_location - the default location to copy files to if there is not a more specific configuration; defaults to empty string
* Boolean verbose - optional - if true, generate verbose log output
* String orchutils_docker_image - required - the name/tag of the mgbpmbiofx/orchutils Docker image
* String gcp_project_id - required - the GCP project id to fetch location connection information from
* String workspace_namespace - optional - the current Terra workspace namespace, used to fetch submission metadata
* String workspace_name - required - the current Terra workspace name, used to fetch location connection information and submission metadata
* String submission_id - optional - the submission id to fetch metadata and logs for

### Output Parameters
* File outputs_manifest - the manifest of outputs and where they were copied to

