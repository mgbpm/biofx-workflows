version 1.0

workflow CopyFileListWorkflow {
    input {
        Array[String] file_list
        Boolean flatten = false
        Boolean recursive = true
        String target_location
        String orchutils_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
        String gcp_project_id
        String workspace_name
    }

    call CopyFilesTask {
        input:
            file_list = file_list,
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


task CopyFilesTask{
    input {
        Array[String] file_list = []
        String? target_location
        Boolean flatten = false
        Boolean recursive = true
        Boolean verbose = false
        String docker_image
        Int disk_size = 75
        String gcp_project_id
        String workspace_name
    }

    command <<<
        set -euxo pipefail
        ROOTDIR="$(pwd)"
        pushd $MGBPMBIOFXPATH/biofx-orchestration-utils

        # if target loc not specified, use ROOTDIR
        FINAL_TARGET_LOC="~{target_location}"
        [ -z "${FINAL_TARGET_LOC}" ] && FINAL_TARGET_LOC="${ROOTDIR}/fetched"

        # run script to setup rclone remotes
        ./bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w "~{workspace_name}" -r "${FINAL_TARGET_LOC}"

        FILE_LIST=(~{sep=" " file_list})
        for src_file in "${FILE_LIST[@]}"
        do
            base_name=$(basename ${src_file})

            # execute script to copy files
            mkdir "${ROOTDIR}/fetched"
            ./bin/copy_files.py ~{if flatten then "--flatten" else ""} ~{if recursive then "" else "--no-recursive"} ~{if verbose then "--verbose" else ""} \
            --source "${src_file}" --target "${FINAL_TARGET_LOC}/${base_name}" \
            --source-files-fofn "${ROOTDIR}/source-file-list.txt" --target-files-fofn "${ROOTDIR}/target-file-list.txt"
        done        

    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        Array[String] source_files = read_lines("source-file-list.txt")
        Array[String] target_files = read_lines("target-file-list.txt")
        Array[File] local_files = glob("fetched/*")
    }
}
