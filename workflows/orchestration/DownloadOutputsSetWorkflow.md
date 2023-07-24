# Download Outputs Workflow
Workflow that distributes workflow output files to various locations as defined by the provided output mapping
configuration.  Intended to be used to download outputs for entity set workflows.  Only the first element of the
default target location list is used; the others are ignored.  The parameter accepts an array to be compatible
with a `this.[setmembers].default_target_location` expression in Terra.

See also `FileUtils.DownloadOutputsTask`.

### Input Parameters
* String outputs_json - the workflow outputs as a stringified JSON object where the keys are the output names
* Array[String] config_json_list - a list of output copy mapping configurations, see `biofx-workflow-configurations` for more details
* String default_target_location_list - default location to copy files to if there is not a more specific configuration; defaults to empty
* Boolean verbose - optional - if true, generate verbose log output
* String orchutils_docker_image - required - the name/tag of the mgbpmbiofx/orchutils Docker image; defaults to "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230719"
* String gcp_project_id - required - the GCP project id to fetch location connection information from; defaults to "mgb-lmm-gcp-infrast-1651079146"
* String workspace_namespace - optional - the current Terra workspace namespace, used to fetch submission metadata
* String workspace_name - required - the current Terra workspace name, used to fetch location connection information and submission metadata
* String submission_id - optional - the submission id to fetch metadata and logs for

### Output Parameters
* File outputs_manifest - the manifest of outputs and where they were copied to

