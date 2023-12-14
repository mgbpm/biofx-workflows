version 1.0

import "../../steps/QCEval.wdl" # for qc of each vcf
import "../../steps/FASTUtils.wdl" # for loading to FAST

workflow SampleLoadingWorkflow {
    input {
        # Dataset and sample inputs
        File input_vcf
        String sample_ID
        String dataset
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker image
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20231129"
        # reference genome
        String reference_build = "GRCh38"
        # qceval inputs
        String qceval_project_type
        String qceval_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20231005"
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

    output {
        # Annotated qceval VCF
        File qceval_vcf_gz = QCEvalTask.output_vcf_gz
        # FAST sample data name
        String fast_sample_data_name = dataset + "_" + sample_ID
        # QCEval load data id
        String qceval_data_load_id = QCEvalLoadTask.data_load_id
    }
}