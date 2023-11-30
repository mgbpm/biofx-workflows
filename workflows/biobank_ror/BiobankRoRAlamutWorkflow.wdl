version 1.0

import "../../steps/VCFUtils.wdl" # for processing input data
import "../../steps/FASTUtils.wdl" # for loading to FAST
import "../../steps/AlamutBatch.wdl" # for annotation
import "../../steps/Utilities.wdl" # for fail task

workflow LoadingAlamutAndGnomadWorkflow {
    input {
        # File for batch annotation
        File batch_annotation_input
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker image
        String orchutils_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20231129"
        # bcftools docker image
        String bcftools_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        # reference genome
        String reference_build = "GRCh38"
        # alamut inputs
        File alamut_db
        File? alamut_fields_tsv
        String alamut_db_name = "alamut_db"
        String alamut_server = "a-ht-na.interactive-biosoftware.com"
        String alamut_port = "80"
        String alamut_user_secret_name = "alamut-batch-ini-user"
        Int alamut_queue_limit = 4
        String alamut_queue_folder = "gs://biofx-task-queue/alamut"
        String alamut_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/alamut:20230630"
        Boolean alamut_save_working_files = false
        String alamut_anno_src_id = "228"
        String alamut_anno_min_age = "P6M"
        # gnomad annotation inputs
        File gnomad_coverage_file
        File gnomad_coverage_file_idx
        Array[String] gnomad_headers
        String gnomad_column_list = "CHROM,POS,INFO/DP_gnomadG"
        # FAST loading inputs
        Boolean has_haploid_sites = false
        String gnomad_data_load_config_name
        String alamut_data_load_config_name = "Alamut"
    }

    # Pre-filter the concat/whole VCF file to input to Alamut to remove
    #  variants that already have Alamut annotations
    call FASTUtils.FASTRemoveAlreadyAnnotatedFromVCFTask {
        input:
            input_vcf = batch_annotation_input,
            output_basename = sub(basename(batch_annotation_input), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".alamutprefilter",
            reference_build = reference_build,
            annotation_source_id = alamut_anno_src_id,
            annotation_min_age = alamut_anno_min_age,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }   
    # Annotate concat/whole VCF with Alamut
    call AlamutBatch.AlamutBatchTask {
        input:
            input_vcf = FASTRemoveAlreadyAnnotatedFromVCFTask.output_vcf_gz,
            output_basename = sub(basename(batch_annotation_input), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".alamut",
            alamut_db = alamut_db,
            alamut_fields_tsv = alamut_fields_tsv,
            alamut_db_name = alamut_db_name,
            alamut_server = alamut_server,
            alamut_port = alamut_port,
            reference_build = reference_build,
            gcp_project_id = gcp_project_id,
            alamut_user_secret_name = alamut_user_secret_name,
            alamut_queue_limit = alamut_queue_limit,
            alamut_queue_folder = alamut_queue_folder,
            docker_image = alamut_docker_image,
            output_working_files = alamut_save_working_files
    }
    # Annotated concat/whole VCF with Gnomad coverage
    call VCFUtils.AnnotateVCFTask as AnnotateGnomadTask {
        input:
            input_vcf = batch_annotation_input,
            output_basename = sub(basename(batch_annotation_input), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".gnomad",
            annotations_file = gnomad_coverage_file,
            annotations_idx_file = gnomad_coverage_file_idx,
            headers_file = write_lines(gnomad_headers),
            column_list = gnomad_column_list,
            docker_image = bcftools_docker_image
    }
    # Load Alamut annotations to FAST
    call FASTUtils.FASTDataLoadTask as AlamutLoadTask {
        input:
            reference_build = reference_build,
            vcf_file = AlamutBatchTask.output_vcf_gz,
            has_haploid_sites = has_haploid_sites,
            data_load_config_name = alamut_data_load_config_name,
            data_load_target = "ANNOTATION_DATA",
            merge_strategy = "MERGE",
            annotation_record_ts = "now",
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    # Load Gnomad annotations to FAST
    call FASTUtils.FASTDataLoadTask as GnomadLoadTask {
        input:
            reference_build = reference_build,
            vcf_file = AnnotateGnomadTask.output_vcf_gz,
            has_haploid_sites = has_haploid_sites,
            data_load_config_name = gnomad_data_load_config_name,
            data_load_target = "ANNOTATION_DATA",
            merge_strategy = "MERGE",
            annotation_record_ts = "now",
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image  
    }
    # Wait for data loads to complete
    call FASTUtils.FASTWaitForDataLoadsTask {
        input:
            data_load_ids = [AlamutLoadTask.data_load_id, GnomadLoadTask.data_load_id],
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
        # Alamut annotated VCF
        File alamut_vcf_gz = AlamutBatchTask.output_vcf_gz
        # gnomAD annotated VCF
        File gnomad_vcf_gz = AnnotateGnomadTask.output_vcf_gz
    }
}