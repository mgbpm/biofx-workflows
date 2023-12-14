version 1.0

import "../../steps/FASTUtils.wdl"
import "../../steps/Utilities.wdl"
import "../../steps/FASTOutputParser.wdl"

workflow FASTParsingWorkflow {
    input {
        # Sample data names in FAST
        String annotated_sample_data_name
        Array[String] fast_data_sample_names
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20231129"
        # reference genome files
        String reference_build = "GRCh38"
        # FAST loading inputs
        Array[String]? fast_annotated_sample_data_regions
        Array[String]? fast_annotated_sample_data_scripts
        String? fast_annotated_sample_data_saved_filter_name
        # Reporting steps
        Boolean create_parsed_output = false
        String fast_parser_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20231206"
        File gil_transcript_exon_count = "gs://lmm-reference-data/annotation/gil_lmm/transcript_exonNum.txt"
        String fast_parser_sample_type = "B"
    }
    
    # Create annotated sample data
    call FASTUtils.FASTCreateAnnotatedSampleDataTask as CreateAnnotatedSampleData {
        input:
            annotated_sample_data_name = annotated_sample_data_name,
            sample_data_names_and_labels = fast_data_sample_names,
            region_names_and_masks = fast_annotated_sample_data_regions,
            scripts = fast_annotated_sample_data_scripts,
            saved_filter_name = fast_annotated_sample_data_saved_filter_name,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    # Wait for annotated sample data to complete
    call FASTUtils.FASTWaitForAnnotatedSampleDataTask as WaitForAnnotatedData {
        input:
            annotated_sample_data_name = CreateAnnotatedSampleData.annotated_sample_data_name_output,
            check_interval_secs = 300,
            max_checks = 72,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    if (WaitForAnnotatedData.wait_result.total != WaitForAnnotatedData.wait_result.completed) {
        call Utilities.FailTask as AnnotatedSampleTimeout {
            input:
                error_message = "FAST annotated sample data creation did not complete within allotted time"
        }
    }
    if (WaitForAnnotatedData.wait_result.total != WaitForAnnotatedData.wait_result.succeeded) {
        call Utilities.FailTask as AnnotatedSampleFailure {
            input:
                error_message = "FAST annotated sample data creation did not succeed"
        }
    }
    if (WaitForAnnotatedData.wait_result.total == WaitForAnnotatedData.wait_result.succeeded) {
        call FASTUtils.FASTExportAnnotatedSampleDataTask {
            input:
                annotated_sample_data_name = CreateAnnotatedSampleData.annotated_sample_data_name_output,
                format = "TXT",
                output_basename = annotated_sample_data_name + ".fastexport",
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }
        if (create_parsed_output) {
            call FASTOutputParser.FASTOutputParserTask {
                input:
                    fast_output_file = FASTExportAnnotatedSampleDataTask.output_file,
                    sample_type = fast_parser_sample_type,
                    reference_build = reference_build,
                    oms_query = "Y",
                    transcript_exonNum = gil_transcript_exon_count,
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name,
                    fast_parser_image = fast_parser_image
            }
        } 
    }

    output {
        # FAST export file
        File? fast_export_file = FASTExportAnnotatedSampleDataTask.output_file
        # FAST Parsed output
        File? fast_parsed_output = FASTOutputParserTask.parsed_report
    }
}