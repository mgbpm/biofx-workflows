version 1.0

import "../../steps/QCEval.wdl"
import "../../steps/FASTUtils.wdl"

workflow BahrainSampleLoadingWorkflow {
    input {
        # Sample data inputs
        File input_vcf
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchutils docker image
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20240625"
        # QCEval inputs
        String qceval_project_type = "WGS"
        String qceval_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20231005"
        # FAST inputs
        String reference_build = "GRCh38"
        Boolean has_haploid_sites = false
        String sample_data_load_config_name = "Sample_VCF_PPM_Eval"
    }

    # Annotate sample VCFs with QCEval INFO field and load annotations into FAST
    call QCEval.QCEvalTask {
        input:
            input_vcf = input_vcf,
            project_type = qceval_project_type,
            docker_image = qceval_docker_image
    }
    call FASTUtils.FASTDataLoadTask as QCEvalLoadTask {
        input:
            reference_build = reference_build,
            vcf_file = QCEvalTask.output_vcf_gz,
            has_haploid_sites = has_haploid_sites,
            sample_data_name = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", ""),
            data_load_config_name = sample_data_load_config_name,
            data_load_target = "SAMPLE_DATA",
            annotation_record_ts = "now",
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    
    output {
        # QCEval VCF
        File qceval_vcf_gz = QCEvalTask.output_vcf_gz
        # FAST sample data name
        String fast_sample_data_name = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
        # QCEval load data id
        String qceval_data_load_id = QCEvalLoadTask.data_load_id
    }
}