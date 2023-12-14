version 1.0

import "../../steps/FASTOutputParser.wdl"

workflow MergeFASTExportsWorkflow {
    input {
        String annotated_sample_data_name
        String sourcename
        # FAST export files (as full paths to bucket locations)
        Array[File] fast_export_files
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # reference genome files
        String reference_build = "GRCh38"
        # Reporting steps
        String python_docker_image = "python:3.10"
        String fast_parser_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20231206"
        File gil_transcript_exon_count = "gs://lmm-reference-data/annotation/gil_lmm/transcript_exonNum.txt"
        String fast_parser_sample_type = "B"
    }
    
    # Merge FAST export files into one large file
    call MergeExportsTask as MergeExports {
        input:
            fast_export_files = fast_export_files,
            docker_image = python_docker_image,
            output_basename = annotated_sample_data_name + ".mergedfastexport",
            sourcename = sourcename
    }
    # Run output parser on merged export output
    call FASTOutputParser.FASTOutputParserTask {
        input:
            fast_output_file = MergeExports.merged_fast_export,
            sample_type = fast_parser_sample_type,
            reference_build = reference_build,
            oms_query = "Y",
            transcript_exonNum = gil_transcript_exon_count,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            fast_parser_image = fast_parser_image
    }

    output {
        # Merged export file
        File merged_export_file = MergeExports.merged_fast_export
        # FAST parsed output
        File? fast_parsed_output = FASTOutputParserTask.parsed_report
    }
}

task MergeExportsTask {
    input {
        Array[File] fast_export_files
        String sourcename
        String output_basename
        String docker_image
        Int mem_size = 2
        Int disk_size = ceil(size(fast_export_files, "GB") * 2.2) + 10
    }

    command <<<
        set -euxo pipefail

        mkdir -p work_dir
        mkdir -p unzipped_export_files

        # Unzip export files
        for file in '~{sep="' '" fast_export_files}'; do
            gunzip -q unzipped_export_files/$file
        done

        # List all unzipped FAST export files in a file
        for file in unzipped_export_files/*; do
             echo ${file} >> work_dir/files_to_merge.txt
        done

        # Merge all unzipped FAST export files
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/merge_export_files.py \
            --files-to-merge work_dir/files_to_merge.txt 
            --sourcename "~{sourcename}" \
            --output-filename "~{output_basename}.txt"
    >>>

    runtime {
        docker: "~{docker_image}"
        memory: "~{mem_size}"
        disks: "local-disk " + disk_size + " SSD"
    }

    output {
        File merged_fast_export = "~{output_basename}.txt"
    }
}