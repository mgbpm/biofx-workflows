version 1.0

task GDHIngestAndFilterTask {
    input {
        String subject_id
        String sample_id
        String gdh_institution = "MGBPM"
        String gdh_project = "Clinical"
        File vcf_file
        String vcf_file_stage_name = "biofx_pipelines"
        String vcf_file_stage_gspath = "gs://gdh-external-stage/biofx_pipelines_nonprod"
        String reference_build = "GRCh38"
        String filter_name_or_code
        String? pipeline_run_id
        Int timeout_minutes = 90
        String gcp_project_id
        String workspace_name
        String docker_image
        Int disk_size = ceil(size(vcf_file, "GB")) + 10
        Boolean verbose = false
    }

    command <<<
        set -euxo pipefail

        # Copy files to GS path that is also GDH external stage
        ./bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w "~{workspace_name}" -r "~{vcf_file_stage_gspath}"
        ./bin/copy_files.py ~{if verbose then "--verbose" else ""} \
            --source "~{vcf_file}" --target "~{vcf_file_stage_gspath}"

        # Start the ingest and filter process
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh \
            -p ~{gcp_project_id} -w ~{workspace_name} -n gdhpipeline > gdhpipeline-client-config.json
        
        cat <<EOF{ 
    "biosample_id1_sys": "LMM-ACCESSION-ID",
    "biosample_id1_id": "~{subject_id}",
    "biosample_id1_lbl": "Accession ID",
    "biosample_id2_sys": "LMM-SAMPLE-ID",
    "biosample_id2_id": "~{sample_id}",
    "biosample_id2_lbl": "Sample ID"
}
EOF > biosample-template.json

        exec_id="~{pipeline_run_id}"
        [ ! -z "$exec_id" ] && exec_id=$(dd if=/dev/random bs=6 count=1 2>>/dev/null | base64 | tr -dC '[:alnum:]')

        echo "~{subject_id}_~{sample_id}-$exec_id" > invoker-execution-id.txt

        $MGBPMBIOFXPATH/biofx-pygdh/bin/run_ingest_and_filter.py ~{if verbose then "--verbose" else ""} \
            --client-config gdhpipeline-client-config.json \
            --run-type "single_sample" \
            --execution-id "~{subject_id}_~{sample_id}-$exec_id" \
            --biosample-default_template "@biosample-template.json" \
            --vcf-file-name "~{vcf_file}" \
            --vcf-file-stage "~{vcf_file_stage_name}" \
            --institution "~{gdh_institution}" \
            --project "~{gdh_project}" \
            --reference-asm "~{reference_build}" \
            --filter-name "~{filter_name_or_code}" \
            --output-file "~{subject_id}_~{sample_id}_~{filter_name_or_code}_matching_variants.json" \
            --timeout-minutes ~{timeout_minutes}
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        String ingest_and_filter_execution_id = read_string("invoker-execution-id.txt")
        File matching_variants = "~{subject_id}_~{sample_id}_~{filter_name_or_code}_matching_variants.json"
    }
}

task GDHOutputParserTask {
    input {
        File gdh_output_file
        String reference_build
        String oms_query = "Y"
        File portable_db_file
        String report_basename = sub(basename(gdh_output_file), "\\.(json|json.gz|JSON|JSON.GZ|txt.gz|txt|TXT.GZ|TXT)$", "")
        String filter_name_or_code
        String gcp_project_id
        String workspace_name
        String gdh_parser_image
    }

    command <<<

        set -uexo pipefail

        # Setup OMS client config
        OMS_PARAM=""
        if [ "~{oms_query}" == "Y" ]; then
            $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh \
                -p ~{gcp_project_id} -w ~{workspace_name} -n oms > oms-client-config.json
            OMS_PARAM="--oms-config oms-client-config.json"
        fi

        # Run the parser
        python gdh-output-parser $OMS_PARAM \
            --input-file "~{gdh_output_file}" \
            --output "~{report_basename}.parsed.txt" \
            --excel-output "~{report_basename}.xlsm" \
            --empty-excel-output "~{report_basename}.xlsx" \
            --gene-info-db "~{portable_db_file}" \
            --build "~{reference_build}" \
            --sample-id-template "{LMM-ACCESSION-ID}_{LMM-SAMPLE-ID}_~{filter_name_or_code}"
    >>>

    runtime {
        docker: "~{gdh_parser_image}"
    }

    output {
        File? parsed_report = report_basename + ".parsed.txt"
        File? nva_report_xlsx = report_basename + ".xlsx"
        File? nva_report_xlsm = report_basename + ".xlsm"
        File nva_report = select_first([nva_report_xlsm, nva_report_xlsx])
    }
}