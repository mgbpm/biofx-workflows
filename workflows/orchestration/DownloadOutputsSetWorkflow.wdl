version 1.0

import "../../steps/FileUtils.wdl"

workflow DownloadOutputsSetWorkflow {
    input {
        String outputs_json
        Array[String] config_json_list
        Array[String] default_target_location_list = []
        String mgbpmbiofx_docker_image
        String gcp_project_id
        String workspace_name
    }

    call FileUtils.DownloadOutputsTask {
        input:
            outputs_json = outputs_json,
            config_json_list = config_json_list,
            default_target_location = if length(default_target_location_list) > 0 then default_target_location_list[0] else "",
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name
    }

    output {
        File outputs_manifest = DownloadOutputsTask.outputs_manifest
    }
}