version 1.0

task FASTDataLoadTask {
    input {
        File vcf_file
        String? default_annotation_src
        String? reference_build
        String data_load_target
        String? sample_data_exists_strategy
        String? merge_strategy
        String? sample_data_coll
        String? sample_data_name
        String? lab_batch_name
        String data_load_config_name
        String? custom_script
        String? annotation_record_ts
        String? email_to
        String gcp_project_id
        String workspace_name
        String docker_image
        Int disk_size = 20
    }

    command <<<
        set -euxo pipefail

        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh \
            -p ~{gcp_project_id} -w ~{workspace_name} -n fast > fast-client-config.json
        
        DEF_ANNO_SRC=""
        [ ! -z "~{default_annotation_src}" ] && DEF_ANNO_SRC="--default-annotation-src ~{default_annotation_src}"

        REF_BUILD=""
        [ ! -z "~{reference_build}" ] && REF_BUILD="--reference-build ~{reference_build}"

        SAMPLE_EXISTS_STRAT=""
        [ ! -z "~{sample_data_exists_strategy}" ] && SAMPLE_EXISTS_STRAT="--sample-data-exists ~{sample_data_exists_strategy}"

        MERGE_STRAT=""
        [ ! -z "~{merge_strategy}" ] && MERGE_STRAT="--merge-strategy ~{merge_strategy}"

        SAMPLE_DATA_COLL=""
        [ ! -z "~{sample_data_coll}" ] && SAMPLE_DATA_COLL="--sample-data-coll ~{sample_data_coll}"

        SAMPLE_DATA_NAME=""
        [ ! -z "~{sample_data_name}" ] && SAMPLE_DATA_NAME="--sample-data-name ~{sample_data_name}"

        LAB_BATCH_NAME=""
        [ ! -z "~{lab_batch_name}" ] && LAB_BATCH_NAME="--lab-batch-name ~{lab_batch_name}"

        CUSTOM_SCRIPT=""
        [ ! -z "~{custom_script}" ] && CUSTOM_SCRIPT="--custom-script ~{custom_script}"

        ANNO_TS=""
        [ ! -z "~{annotation_record_ts}" -a "now" != "~{annotation_record_ts}" ] && ANNO_TS="--annotation-record-ts ~{annotation_record_ts}"
        [ "now" == "~{annotation_record_ts}" ] && ANNO_TS="--annotation-record-ts $(date -Iseconds)"

        EMAIL_TO=""
        [ ! -z "~{email_to}" ] && EMAIL_TO="--email-to \"~{email_to}\""

        OPTIONAL_PARAMS="$DEF_ANNO_SRC $REF_BUILD $SAMPLE_EXISTS_STRAT $MERGE_STRAT"
        OPTIONAL_PARAMS="$OPTIONAL_PARAMS $SAMPLE_DATA_COLL $SAMPLE_DATA_NAME $LAB_BATCH_NAME $CUSTOM_SCRIPT $ANNO_TS $EMAIL_TO"

        $MGBPMBIOFXPATH/biofx-pyfast/bin/upload_data.py \
            --client-config fast-client-config.json \
            --vcf-file "~{vcf_file}" --target ~{data_load_target} \
            --load-config-name "~{data_load_config_name}" $OPTIONAL_PARAMS  > data-load-id.txt
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        String data_load_id = read_string("data-load-id.txt")
    }
}

struct FASTWaitResult {
    Int total
    Int completed
    Int succeeded
    Int failed
    Int canceled
}

task FASTWaitForDataLoadsTask {
    input {
        Array[String] data_load_ids
        Int check_interval_secs = 60
        Int max_checks = 24 * 60
        String gcp_project_id
        String workspace_name
        String docker_image
    }

    command <<<
        set -euxo pipefail

        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh \
            -p ~{gcp_project_id} -w ~{workspace_name} -n fast > fast-client-config.json

        $MGBPMBIOFXPATH/biofx-pyfast/bin/wait_for_data_loads.py \
            --client-config fast-client-config.json \
            --data-load-ids ~{sep="," data_load_ids} \
            --check-interval-secs ~{check_interval_secs} --max-checks ~{max_checks} | tee wait-result.txt
    >>>

    runtime {
        docker: "~{docker_image}"
    }

    output {
        # The linter complains about this line, but it works as expected
        FASTWaitResult wait_result = read_json("wait-result.txt")
    }
}

task FASTCreateAnnotatedSampleDataTask {
    input {
        String annotated_sample_data_name
        Array[String] sample_data_names_and_labels
        Array[String]? region_names_and_masks
        Array[String]? scripts
        String? saved_filter_name
        String gcp_project_id
        String workspace_name
        String docker_image
        Int disk_size = 20
    }

    command <<<
        set -euxo pipefail

        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh \
            -p ~{gcp_project_id} -w ~{workspace_name} -n fast > fast-client-config.json

        REGIONS=""
        [ ! -z "~{sep='' region_names_and_masks}" ] && REGIONS="--regions ~{sep=',' region_names_and_masks}"

        SCRIPTS=""
        [ ! -z "~{sep='' scripts}" ] && SCRIPTS="--scripts ~{sep=',' scripts}"

        SAVEDFILTER=""
        [ ! -z "~{saved_filter_name}" ] && SAVEDFILTER="--saved-filter ~{saved_filter_name}"

        $MGBPMBIOFXPATH/biofx-pyfast/bin/create_annotated_sample_data.py \
            --client-config fast-client-config.json \
            --name "~{annotated_sample_data_name}" \
            --sample-datas "~{sep=',' sample_data_names_and_labels}" \
            $REGIONS $SCRIPTS $SAVEDFILTER

    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        String annotated_sample_data_name_output = annotated_sample_data_name
    }
}

task FASTWaitForAnnotatedSampleDataTask {
    input {
        String annotated_sample_data_name
        Int check_interval_secs = 60
        Int max_checks = 24 * 60
        String gcp_project_id
        String workspace_name
        String docker_image
    }

    command <<<
        set -euxo pipefail
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh \
            -p ~{gcp_project_id} -w ~{workspace_name} -n fast > fast-client-config.json

        $MGBPMBIOFXPATH/biofx-pyfast/bin/wait_for_annotated_sample_data.py \
            --client-config fast-client-config.json \
            --name ~{annotated_sample_data_name} \
            --check-interval-secs ~{check_interval_secs} --max-checks ~{max_checks} > wait-result.txt

        WAIT_RESULT=$(cat wait-result.txt | tr [A-Z] [a-z])

    >>>

    runtime {
        docker: "~{docker_image}"
    }

    output {
        # The linter complains about this line, but it works as expected
        FASTWaitResult wait_result = read_json("wait-result.txt")
    }
}

task FASTRemoveAlreadyAnnotatedFromVCFTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".alreadyann"
        String reference_build
        String? annotation_source_id
        String? annotation_field_id
        String? annotation_min_age
        String gcp_project_id
        String workspace_name
        String docker_image
        Int batch_size = 20
    }

    command <<<
        set -euxo pipefail
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh \
            -p ~{gcp_project_id} -w ~{workspace_name} -n fast > fast-client-config.json

        ANNOTATION_SRC_PARAM=""
        [ ! -z "~{annotation_source_id}" ] && ANNOTATION_SRC_PARAM="--annotation-src-id ~{annotation_source_id}"

        ANNOTATION_FIELD_PARAM=""
        [ ! -z "~{annotation_field_id}" ] && ANNOTATION_SRC_PARAM="--annotation-field-id ~{annotation_field_id}"

        ANNOTATION_MIN_AGE_PARAM=""
        [ ! -z "~{annotation_min_age}" ] && ANNOTATION_SRC_PARAM="--annotation-min-age ~{annotation_min_age}"

        $MGBPMBIOFXPATH/biofx-pyfast/bin/filter_already_annotated.py \
            --client-config fast-client-config.json \
            --vcf-file "~{input_vcf}" --vcf-output-file "~{output_basename}.vcf.gz" \
            --reference-build ~{reference_build} --batch-size ~{batch_size} \
            $ANNOTATION_SRC_PARAM $ANNOTATION_FIELD_PARAM $ANNOTATION_MIN_AGE_PARAM
    >>>

    runtime {
        docker: "~{docker_image}"
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
    }
}