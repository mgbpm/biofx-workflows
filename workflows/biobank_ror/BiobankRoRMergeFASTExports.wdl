version 1.0

import "../../steps/FASTOutputParser.wdl"

workflow MergeFASTExportsWorkflow {
    input {
        # Sample data name for merged exports
        String annotated_sample_data_name
        # FAST export files (as full paths to bucket locations)
        Array[File] fast_export_files
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Reference genome files
        String reference_build = "GRCh38"
        # Reporting steps
        String fast_parser_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20250111"
        File portable_db_file = "gs://lmm-reference-data/annotation/gil_lmm/gene_info.db"
        String fast_parser_sample_type = "B"
        Boolean gatk_source = false
    }
    
    # Merge FAST export files into one large file
    call MergeExportsTask as MergeExports {
        input:
            fast_export_files = fast_export_files,
            docker_image = fast_parser_image,
            output_basename = annotated_sample_data_name + ".mergedfastexport",
            sourcename = annotated_sample_data_name
    }
    # Run output parser on merged export output
    call FASTOutputParser.FASTOutputParserTask {
        input:
            fast_output_file = MergeExports.merged_fast_export,
            sample_type = fast_parser_sample_type,
            reference_build = reference_build,
            oms_query = "Y",
            portable_db_file = portable_db_file,
            gatk_source = gatk_source,
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