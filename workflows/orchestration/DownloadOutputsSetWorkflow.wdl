version 1.0

import "../../steps/FileUtils.wdl"

workflow DownloadOutputsSetWorkflow {
    input {
        String outputs_json
        Array[String?] config_json_list
        String? default_target_location
        Boolean verbose = false
        String orchutils_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String? workspace_namespace
        String workspace_name
        String? submission_id
    }

    call FileUtils.DownloadOutputsTask {
        input:
            outputs_json = outputs_json,
            config_json_list = select_all(config_json_list),
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