version 1.0

struct FileMatcher {
    String file_type
    Array[String]? path_matches
    Boolean? required
}

task CopyFilesTask {
    input {
        String source_location
        Array[String] file_types = []
        Array[String] file_match_keys = []
        Array[FileMatcher] file_matchers = []
        String? target_location
        Boolean flatten = false
        Boolean recursive = true
        Boolean verbose = false
        String docker_image
        Int disk_size = 75
        String gcp_project_id
        String workspace_name
    }

    Boolean has_file_matchers = defined(file_matchers) && length(file_matchers) > 0

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
        [ -z "${FILE_TYPE_LIST}" -a "~{has_file_matchers}" == "false" ] && FILE_TYPE_ARG="--optional-file-types-all"
        [ ! -z "${FILE_TYPE_LIST}" ] && FILE_TYPE_ARG="--optional-file-types ${FILE_TYPE_LIST}"
        MATCH_KEY_LIST="~{sep="," file_match_keys}"
        [ ! -z "${MATCH_KEY_LIST}" ] && MATCH_KEY_ARG="--filter-keys ${MATCH_KEY_LIST}"

        FILE_MATCHERS_ARG=""
        [ "~{has_file_matchers}" == "true" ] && FILE_MATCHERS_ARG="--file-matchers ~{write_json(file_matchers)}"

        # execute script to copy files
        mkdir "${ROOTDIR}/fetched"
        ./bin/copy_files.py ~{if flatten then "--flatten" else ""} ~{if recursive then "" else "--no-recursive"} ~{if verbose then "--verbose" else ""} \
            --source "~{source_location}" --target "${FINAL_TARGET_LOC}" ${FILE_TYPE_ARG} ${MATCH_KEY_ARG} ${FILE_MATCHERS_ARG} \
            --source-files-fofn "${ROOTDIR}/source-file-list.txt" --target-files-fofn "${ROOTDIR}/target-file-list.txt"

        # test if target files are local
        touch "${ROOTDIR}/local-file-list.txt"
        FIRST_TARGET_PATH=$(head -1 "${ROOTDIR}/target-file-list.txt")
        [ ! -z "${FIRST_TARGET_PATH}" -a -f "${FIRST_TARGET_PATH}" ] && cp "${ROOTDIR}/target-file-list.txt" "${ROOTDIR}/local-file-list.txt"
            
        popd
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

task FetchFilesTask {
    input {
        String data_location
        Boolean recursive = true
        Array[String] file_types = []
        Array[String] file_match_keys = []
        Array[FileMatcher] file_matchers = []
        Boolean verbose = false
        String docker_image
        Int disk_size = 75
        String gcp_project_id
        String workspace_name
        File? empty_output_placeholder
    }

    Boolean has_file_matchers = defined(file_matchers) && length(file_matchers) > 0

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
        [ -z "${FILE_TYPE_LIST}" -a "~{has_file_matchers}" == "false" ] && FILE_TYPE_ARG="--optional-file-types-all"
        [ ! -z "${FILE_TYPE_LIST}" ] && FILE_TYPE_ARG="--optional-file-types ${FILE_TYPE_LIST}"
        MATCH_KEY_LIST="~{sep="," file_match_keys}"
        [ ! -z "${MATCH_KEY_LIST}" ] && MATCH_KEY_ARG="--filter-keys ${MATCH_KEY_LIST}"

        FILE_MATCHERS_ARG=""
        [ "~{has_file_matchers}" == "true" ] && FILE_MATCHERS_ARG="--file-matchers ~{write_json(file_matchers)}"

        # execute script to copy files
        mkdir "${ROOTDIR}/fetched"
        ./bin/copy_files.py \
            --source "~{data_location}" --target "${ROOTDIR}/fetched" --flatten ${FILE_TYPE_ARG} ${MATCH_KEY_ARG} ${FILE_MATCHERS_ARG} ~{if verbose then "--verbose" else ""} \
            --source-files-fofn "${ROOTDIR}/source-file-list.txt" --target-files-fofn "${ROOTDIR}/target-file-list.txt" \
            ~{if recursive then "" else "--no-recursive"}
        popd

        # extract specific file types from list for outputs
        set +e
        grep -i "[.]bam$" target-file-list.txt | sort | head -1 | xargs basename > target-file-list-bam.txt
        grep -i "[.]bai$" target-file-list.txt | sort | head -1 | xargs basename > target-file-list-bai.txt
        grep -i "[.]cram$" target-file-list.txt | sort | head -1 | xargs basename > target-file-list-cram.txt
        grep -i "[.]crai$" target-file-list.txt | sort | head -1 | xargs basename > target-file-list-crai.txt
        grep -iE "[.](vcf|vcf.gz|vcf.bgz|vcf.bz2|gvcf|gvcf.gz|gvcf.bgz|gvcf.bz2)$" target-file-list.txt | sort | head -1 | xargs basename > target-file-list-vcf.txt
        grep -iE "[.](vcf|gvcf)" target-file-list.txt | grep -iE "[.](tbi|idx|csi)$" | sort | head -1 | xargs basename > target-file-list-vcfidx.txt
        grep -iE "[.](bcf|bcf.gz)$" target-file-list.txt | sort | head -1 | xargs basename > target-file-list-bcf.txt

        # diagnostic output for debugging
        ls -al "${ROOTDIR}/fetched"
        cat target-file-list*
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        Array[File] all_files = glob("fetched/*")
        String glob_dir = sub(if length(all_files) > 0 then all_files[0] else "", "[^/]*$", "")
        File? bam = if size("target-file-list-bam.txt") > 0 then glob_dir + read_string("target-file-list-bam.txt") else empty_output_placeholder
        File? bai = if size("target-file-list-bai.txt") > 0 then glob_dir + read_string("target-file-list-bai.txt") else empty_output_placeholder
        File? cram = if size("target-file-list-cram.txt") > 0 then glob_dir + read_string("target-file-list-cram.txt") else empty_output_placeholder
        File? crai = if size("target-file-list-crai.txt") > 0 then glob_dir + read_string("target-file-list-crai.txt") else empty_output_placeholder
        File? vcf = if size("target-file-list-vcf.txt") > 0 then glob_dir + read_string("target-file-list-vcf.txt") else empty_output_placeholder
        File? vcf_index = if size("target-file-list-vcfidx.txt") > 0 then glob_dir + read_string("target-file-list-vcfidx.txt") else empty_output_placeholder
        File? bcf = if size("target-file-list-bcf.txt") > 0 then glob_dir + read_string("target-file-list-bcf.txt") else empty_output_placeholder
    }
}

task DownloadOutputsTask {
    input {
        String outputs_json
        Array[String] config_json_list
        String? default_target_location
        Boolean verbose = false
        String docker_image
        String gcp_project_id
        String? workspace_namespace
        String workspace_name
        String? submission_id
    }

    command <<<
        set -euxo pipefail
        ROOTDIR="$(pwd)"
        pushd $MGBPMBIOFXPATH/biofx-orchestration-utils

        DEF_TARGET_ARG=""
        [ ! -z "~{default_target_location}" ] && DEF_TARGET_ARG="--default-target-location ~{default_target_location}"

        WSNS_ARG=""
        [ ! -z "~{workspace_namespace}" ] && WSNS_ARG="--workspace-namespace ~{workspace_namespace}"

        SUBID_ARG=""
        [ ! -z "~{submission_id}" ] && SUBID_ARG="--submission-id ~{submission_id}"

        # run script to setup rclone remotes
        REMOTES=$(./bin/get_outputs_remotes.py --outputs "~{write_lines([outputs_json])}" --config "~{write_lines(config_json_list)}" ${DEF_TARGET_ARG})
        for remote in ${REMOTES}
        do
            ./bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w "~{workspace_name}" -n ${remote}
        done

        # execute script to copy files
        ./bin/copy_outputs.py --outputs "~{write_lines([outputs_json])}" --config "~{write_lines(config_json_list)}" ~{if verbose then "--verbose" else ""} \
            ${DEF_TARGET_ARG} ${WSNS_ARG} --workspace-name "~{workspace_name}" ${SUBID_ARG} --local-manifest-file "${ROOTDIR}/copy-manifest.json"
        popd
    >>>

    runtime {
        docker: "~{docker_image}"
    }

    output {
        File outputs_manifest = "copy-manifest.json"
    }
}
