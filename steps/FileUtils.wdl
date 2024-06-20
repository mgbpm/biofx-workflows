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

task GCPCopyAndRenameVCF {
    input {
        String source_file
        String sample_id
        String subject_id
        String target_location
        String docker_image
        Int disk_size = 75
    }

    command <<<
        set -euxo pipefail
        #run gsutil cp to move file from one bucket to another
        gsutil -q stat "~{source_file}"
        PATH_EXIST=$?
        if [ ${PATH_EXIST} -eq 0 ]; then
            echo "~{source_file} exists."
        else
            echo "~{source_file} does not exist."
            exit 1
        fi

        #copy file from one bucket to another
        FILE_NAME=$(basename "~{source_file}")
        echo "FILE_NAME variable: ${FILE_NAME}"
        #strip off ending 
        #subject_id="${FILE_NAME%.*}"
        echo "subject_id variable: ~{subject_id}"
        #grab extension
        EXTENSION=$(echo "$FILE_NAME" | sed 's/^[^.]*\.//')
        echo "EXTENSION variable: ${EXTENSION}"
        #set new filename
        SAMPLE_ID_CLEAN=$(sed -e 's/^"//' -e 's/"$//' <<<"~{sample_id}")
        FILE_NAME_NEW="~{subject_id}_$SAMPLE_ID_CLEAN.$EXTENSION"
        echo "renamed ${FILE_NAME} to ${FILE_NAME_NEW}"
        gsutil cp "~{source_file}" "~{target_location}"/${FILE_NAME_NEW}
        COPY_STATUS_VCF=$?
        #download vcf locally
        gsutil cp "~{source_file}" .
        #generate index file
        tabix -p vcf ${FILE_NAME_NEW}
        #copy index file
        gsutil cp ${FILE_NAME_NEW}.tbi "~{target_location}/${FILE_NAME_NEW}.tbi"
        COPY_STATUS_TBI=$?
        #relocalize renamed vcf
        gsutil cp "~{target_location}"/${FILE_NAME_NEW} .

        #check copy success
        if [ ${COPY_STATUS_VCF} -eq 0 ] && [ ${COPY_STATUS_TBI} -eq 0 ]; then
            echo 'Successfully copied "~{source_file}" to "~{target_location}"/${FILE_NAME_NEW}.'
            echo 'Successfully copied "~{source_file}" to "~{target_location}"/${FILE_NAME_NEW}.' > copy-manifest.log
        else
            echo "Unsuccessfull copy. Error code ${COPY_STATUS}"
            echo "Unsuccessfull copy. Error code ${COPY_STATUS}" > copy-manifest.log
            exit ${COPY_STATUS_VCF}
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File outputs_manifest = "copy-manifest.log"
        File output_vcf = "~{subject_id}_~{sample_id}.vcf.gz"
        File output_vcf_index = "~{subject_id}_~{sample_id}.vcf.gz.tbi"
    }
}
