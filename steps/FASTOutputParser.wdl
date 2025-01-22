version 1.0

task FASTOutputParserTask {
    input {
        File fast_output_file
        String sample_type
        String reference_build = "GRCh38"
        String oms_query = "Y"
        File portable_db_file
        String report_basename = sub(basename(fast_output_file), "\\.(txt.gz|txt|TXT.GZ|TXT)$", "")
        Boolean gatk_source = false
        String gcp_project_id
        String workspace_name
        String fast_parser_image
        Int preemptible = 1
    }

    Int cpu = 1
    String file_name = basename(fast_output_file)
    String file_basename = sub(file_name, "\\.(txt.gz|txt|TXT.GZ|TXT)$", "")

    runtime {
        docker: "~{fast_parser_image}"
        cpu: "~{cpu}"
        memory: "1 GB"
        preemptible: preemptible
    }

    command <<<
        set -uexo pipefail

        # Setup OMS client config
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh \
            -p ~{gcp_project_id} -w ~{workspace_name} -n oms > oms-client-config.json

        # Move input file to working dir so that output files are predictably located
        mv "~{fast_output_file}" "~{file_name}"

        # Run the parser
        if [ "~{gatk_source}" == "true" ]; then
            $MGBPMBIOFXPATH/biofx-fast-output-parser/bin/run_parser.py -f "~{file_name}" \
                        -s "~{sample_type}" \
                        -o "~{oms_query}" \
                        -e "~{portable_db_file}" \
                        -b "~{reference_build}" \
                        -k oms-client-config.json \
                        -a
        else
           $MGBPMBIOFXPATH/biofx-fast-output-parser/bin/run_parser.py -f "~{file_name}" \
                        -s "~{sample_type}" \
                        -o "~{oms_query}" \
                        -e "~{portable_db_file}" \
                        -b "~{reference_build}" \
                        -k oms-client-config.json 
        fi

        # Rename the output file to match the report basename
        if [[ "~{file_basename}" != "~{report_basename}" ]]
        then
            if [ -f "~{file_basename}.xlsx" ]
            then
                mv "~{file_basename}.xlsx" "~{report_basename}.xlsx"
            fi
            if [ -f "~{file_basename}.xlsm" ]
            then
                mv "~{file_basename}.xlsm" "~{report_basename}.xlsm"
            fi
        fi
    >>>

    output {
        File? parsed_report = file_basename + ".parsed.txt"
        File? nva_report_xlsx = report_basename + ".xlsx"
        File? nva_report_xlsm = report_basename + ".xlsm"
        File nva_report = select_first([nva_report_xlsm, nva_report_xlsx])
    }
}

task MergeFASTExportsTask {
    input {
        Array[File] fast_export_files
        String sourcename
        String output_basename
        String docker_image
        Int disk_size = ceil(size(fast_export_files, "GB") * 2.2) + 10
    }

    command <<<
        set -euxo pipefail

        # List all unzipped FAST export files in a file
        for c in '~{sep="' '" fast_export_files}'; do
            echo $c >> merge_these_files.txt
        done

        # Merge all FAST export files
        $MGBPMBIOFXPATH/biofx-fast-output-parser/bin/merge_export_files.py \
            --files-to-merge merge_these_files.txt --sourcename "~{sourcename}" --output-filename "~{output_basename}.txt"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
    }

    output {
        File merged_fast_export = "~{output_basename}.txt"
    }
}