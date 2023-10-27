version 1.0

import "../../steps/QCEval.wdl" # for qc of each vcf
import "../../steps/FASTUtils.wdl" # for loading to FAST
import "../../steps/Utilities.wdl" # for fail task

workflow SampleLoadingWorkflow {
    input {
        # Input file (result from data prep workflow)
        File input_vcf
        String sample_ID
        # About the dataset
        String dataset
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker image
        String orchutils_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230921"
        # reference genome
        String reference_build = "GRCh38"
        # qceval inputs
        String qceval_project_type
        String qceval_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20231005"
        # FAST loading inputs
        Boolean has_haploid_sites = false
        String sample_data_load_config_name = "Sample_VCF_PPM_Eval"
    }

    # Annotate individual sample VCF with QCEval INFO field
    call QCEval.QCEvalTask {
        input:
            input_vcf = input_vcf,
            project_type = qceval_project_type,
            output_basename = dataset + "_" + sample_ID + ".qceval",
            docker_image = qceval_docker_image
    }
    # Load sample data to FAST
    call FASTUtils.FASTDataLoadTask as QCEvalLoadTask {
        input:
            reference_build = reference_build,
            vcf_file = QCEvalTask.output_vcf_gz,
            has_haploid_sites = has_haploid_sites,
            sample_data_name = dataset + "_" + sample_ID,
            data_load_config_name = sample_data_load_config_name,
            data_load_target = "SAMPLE_DATA",
            annotation_record_ts = "now",
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    # Wait for data loads to complete
    call FASTUtils.FASTWaitForDataLoadsTask {
        input:
            data_load_ids = [QCEvalLoadTask.data_load_id],
            check_interval_secs = 300,
            max_checks = 72,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    if (FASTWaitForDataLoadsTask.wait_result.total != FASTWaitForDataLoadsTask.wait_result.completed) {
        call Utilities.FailTask as DataLoadTimeout {
            input:
                error_message = "One or more FAST data loads did not complete within allotted time"
        }
    }
    if (FASTWaitForDataLoadsTask.wait_result.total != FASTWaitForDataLoadsTask.wait_result.succeeded) {
        call Utilities.FailTask as DataLoadFailure {
            input:
                error_message = "One or more FAST data loads did not succeed"
        }
    }

    output {
        # Annotated qceval VCF
        File qceval_vcf_gz = QCEvalTask.output_vcf_gz
        # FAST sample data name
        String fast_sample_data_name = dataset + "_" + sample_ID
    }
}