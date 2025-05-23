version 1.0

import "../../steps/FileUtils.wdl"

workflow DownloadOutputsWorkflow {
    input {
        String outputs_json
        String? config_json
        String? default_target_location
        Boolean verbose = false
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String? workspace_namespace
        String workspace_name
        String? submission_id
    }

    call FileUtils.DownloadOutputsTask {
        input:
            outputs_json = outputs_json,
            config_json_list = if defined(config_json) then [select_first([config_json])] else [],
            default_target_location = default_target_location,
            verbose = verbose,
            docker_image = orchutils_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_namespace = workspace_namespace,
            workspace_name = workspace_name,
            submission_id = submission_id
    }

    output {
        File outputs_manifest = DownloadOutputsTask.outputs_manifest
    }
}