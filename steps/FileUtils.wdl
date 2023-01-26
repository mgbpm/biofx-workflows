version 1.0

task CopyFilesTask {
    input {
        String source_location
        Array[String] file_types = []
        Array[String] file_match_keys = []
        String? target_location
        Boolean flatten = false
        Boolean recursive = true
        String mgbpmbiofx_docker_image
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
        ./bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w "~{workspace_name}" -r "~{source_location}"
        ./bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w "~{workspace_name}" -r "${FINAL_TARGET_LOC}"

        # build inputs for copy files script
        FILE_TYPE_ARG=""
        MATCH_KEY_ARG=""
        FILE_TYPE_LIST="~{sep="," file_types}"
        [ -z "${FILE_TYPE_LIST}" ] && FILE_TYPE_ARG="--optional-file-types-all"
        [ ! -z "${FILE_TYPE_LIST}" ] && FILE_TYPE_ARG="--optional-file-types ${FILE_TYPE_LIST}"
        MATCH_KEY_LIST="~{sep="," file_match_keys}"
        [ ! -z "${MATCH_KEY_LIST}" ] && MATCH_KEY_ARG="--filter-keys ${MATCH_KEY_LIST}"

        # execute script to copy files
        mkdir "${ROOTDIR}/fetched"
        ./bin/copy_files.py ~{if flatten then "--flatten" else ""} ~{if recursive then "" else "--no-recursive"} --verbose \
            --source "~{source_location}" --target "${FINAL_TARGET_LOC}" ${FILE_TYPE_ARG} ${MATCH_KEY_ARG} \
            --source-files-fofn "${ROOTDIR}/source-file-list.txt" --target-files-fofn "${ROOTDIR}/target-file-list.txt"

        # test if target files are local
        touch "${ROOTDIR}/local-file-list.txt"
        FIRST_TARGET_PATH=$(head -1 "${ROOTDIR}/target-file-list.txt")
        [ ! -z "${FIRST_TARGET_PATH}" -a -f "${FIRST_TARGET_PATH}" ] && cp "${ROOTDIR}/target-file-list.txt" "${ROOTDIR}/local-file-list.txt"
            
        popd
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "4GB"
    }

    output {
        Array[String] source_files = read_lines("source-file-list.txt")
        Array[String] target_files = read_lines("target-file-list.txt")
        Array[File] local_files = glob("fetched/*")
    }
}

task FetchFilesTask {
    input {
        String data_location
        Boolean recursive = true
        Array[String] file_types = []
        Array[String] file_match_keys = []
        String mgbpmbiofx_docker_image
        String gcp_project_id
        String workspace_name
        File? empty_output_placeholder
    }

    command <<<
        set -euxo pipefail
        ROOTDIR="$(pwd)"
        pushd $MGBPMBIOFXPATH/biofx-orchestration-utils

        # run script to setup rclone remotes
        ./bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w "~{workspace_name}" -r "~{data_location}"

        # build inputs for copy files script
        FILE_TYPE_ARG=""
        MATCH_KEY_ARG=""
        FILE_TYPE_LIST="~{sep="," file_types}"
        [ -z "${FILE_TYPE_LIST}" ] && FILE_TYPE_ARG="--optional-file-types-all"
        [ ! -z "${FILE_TYPE_LIST}" ] && FILE_TYPE_ARG="--optional-file-types ${FILE_TYPE_LIST}"
        MATCH_KEY_LIST="~{sep="," file_match_keys}"
        [ ! -z "${MATCH_KEY_LIST}" ] && MATCH_KEY_ARG="--filter-keys ${MATCH_KEY_LIST}"

        # execute script to copy files
        mkdir "${ROOTDIR}/fetched"
        ./bin/copy_files.py \
            --source "~{data_location}" --target "${ROOTDIR}/fetched" --flatten ${FILE_TYPE_ARG} ${MATCH_KEY_ARG} --verbose \
            --source-files-fofn "${ROOTDIR}/source-file-list.txt" --target-files-fofn "${ROOTDIR}/target-file-list.txt" \
            ~{if recursive then "" else "--no-recursive"}
        popd

        # extract specific file types from list for outputs
        set +e
        grep -i "[.]bam$" target-file-list.txt > target-file-list-bam.txt
        grep -i "[.]bai$" target-file-list.txt > target-file-list-bai.txt
        grep -i "[.]cram$" target-file-list.txt > target-file-list-cram.txt
        grep -i "[.]crai$" target-file-list.txt > target-file-list-crai.txt
        grep -iE "[.](vcf|vcf.gz|vcf.bgz|vcf.bz2|gvcf|gvcf.gz|gvcf.bgz|gvcf.bz2)$" target-file-list.txt > target-file-list-vcf.txt
        grep -iE "[.](vcf|gvcf)" target-file-list.txt | grep -iE "[.](tbi|idx)$" > target-file-list-vcfidx.txt

        # diagnostic output for debugging
        ls -al "${ROOTDIR}/fetched"
        cat target-file-list*
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "4GB"
    }

    output {
        File? bam = if size("target-file-list-bam.txt") > 0 then read_string("target-file-list-bam.txt") else empty_output_placeholder
        File? bai = if size("target-file-list-bai.txt") > 0 then read_string("target-file-list-bai.txt") else empty_output_placeholder
        File? cram = if size("target-file-list-cram.txt") > 0 then read_string("target-file-list-cram.txt") else empty_output_placeholder
        File? crai = if size("target-file-list-crai.txt") > 0 then read_string("target-file-list-crai.txt") else empty_output_placeholder
        File? vcf = if size("target-file-list-vcf.txt") > 0 then read_string("target-file-list-vcf.txt") else empty_output_placeholder
        File? vcf_index = if size("target-file-list-vcfidx.txt") > 0 then read_string("target-file-list-vcfidx.txt") else empty_output_placeholder
        Array[File] all_files = glob("fetched/*")
    }
}

task DownloadOutputsTask {
    input {
        String outputs_json
        String config_json
        String? default_target_location
        String mgbpmbiofx_docker_image
        String gcp_project_id
        String workspace_name
    }

    command <<<
        set -euxo pipefail
        ROOTDIR="$(pwd)"
        pushd $MGBPMBIOFXPATH/biofx-orchestration-utils

        DEF_TARGET_ARG=""
        [ ! -z "~{default_target_location}" ] && DEF_TARGET_ARG="--default-target-location ~{default_target_location}"

        # run script to setup rclone remotes
        REMOTES=$(./bin/get_outputs_remotes.py --outputs "~{write_lines([outputs_json])}" --config "~{write_lines([config_json])}" ${DEF_TARGET_ARG})
        for remote in ${REMOTES}
        do
            ./bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w "~{workspace_name}" -n ${remote}
        done

        # execute script to copy files
        ./bin/copy_outputs.py --outputs "~{write_lines([outputs_json])}" --config "~{write_lines([config_json])}" --verbose \
            ${DEF_TARGET_ARG} --local-manifest-file "${ROOTDIR}/copy-manifest.json"
        popd
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "4GB"
    }

    output {
        File outputs_manifest = "copy-manifest.json"
    }
}
