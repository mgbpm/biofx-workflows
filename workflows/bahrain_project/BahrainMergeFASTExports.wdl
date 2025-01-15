version 1.0

import "../../steps/FASTOutputParser.wdl"

workflow BahrainMergeExportsWorkflow {
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
        String fast_parser_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20240130"
        File gil_transcript_exon_count = "gs://lmm-reference-data/annotation/gil_lmm/transcript_exonNum.txt"
        String fast_parser_sample_type
        Boolean gatk_source = false
    }
    
    # Merge FAST export files into one large file
    call FASTOutputParser.MergeFASTExportsTask {
        input:
            fast_export_files = fast_export_files,
            docker_image = fast_parser_image,
            output_basename = annotated_sample_data_name + ".mergedfastexport",
            sourcename = annotated_sample_data_name
    }
    # Run output parser on merged export output
    call FASTOutputParser.FASTOutputParserTask {
        input:
            fast_output_file = MergeFASTExportsTask.merged_fast_export,
            sample_type = fast_parser_sample_type,
            reference_build = reference_build,
            oms_query = "Y",
            transcript_exonNum = gil_transcript_exon_count,
            gatk_source = gatk_source,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            fast_parser_image = fast_parser_image
    }

    output {
        # Merged export file
        File merged_export_file = MergeFASTExportsTask.merged_fast_export
        # FAST parsed output
        File? fast_parsed_output = FASTOutputParserTask.parsed_report
    }
}