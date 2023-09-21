version 1.0

import "../../steps/FileUtils.wdl"

workflow CopySampleFilesWorkflow {
    input {
        String sample_id
        String source_location
        Boolean flatten = false
        Boolean recursive = true
        String target_location
        String orchutils_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
        String gcp_project_id
        String workspace_name
    }

    call FileUtils.CopyFilesTask {
        input:
            source_location = source_location,
            file_match_keys = [sample_id],
            target_location = target_location,
            flatten = flatten,
            recursive = recursive,
            docker_image = orchutils_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name
    }

    output {
        Array[String] source_files = CopyFilesTask.source_files
        Array[String] target_files = CopyFilesTask.target_files
        Array[File] local_files = CopyFilesTask.local_files
        String new_sample_data_location = target_location
    }
}