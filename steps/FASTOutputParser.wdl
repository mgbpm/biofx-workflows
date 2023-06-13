version 1.0

task FASTOutputParserTask {
    input {
        File fast_output_file
        String sample_type
        String reference_build = "GRCh38"
        String oms_query = "Y"
        File transcript_exonNum
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
        $MGBPMBIOFXPATH/biofx-fast-output-parser/bin/run_parser.py -f "~{file_name}" \
                    -s "~{sample_type}" \
                    -o "~{oms_query}" \
                    -e "~{transcript_exonNum}" \
                    -b "~{reference_build}" \
                    -k oms-client-config.json
    >>>

    output {
        File? parsed_report = file_basename + ".parsed.txt"
        File? nva_report_xlsx = file_basename + ".xlsx"
        File? nva_report_xlsm = file_basename + ".xlsm"
        File nva_report = select_first([nva_report_xlsm, nva_report_xlsx])
    }
}