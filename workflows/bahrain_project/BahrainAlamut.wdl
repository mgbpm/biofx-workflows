version 1.0

import "../../steps/AlamutBatch.wdl"
import "../../steps/FASTUtils.wdl"

workflow BahrainAlamutWorkflow {
    input {
        # Collective or merged vcf
		File merged_vcf
		# Batch name
		String batch_name
		# Reference genome files
        String reference_build = "GRCh38"
		# GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
		# Orchutils docker image
		String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20240625"
		# Alamut inputs
        File alamut_db
        File? alamut_fields_tsv
        String alamut_db_name = "alamut_db"
        String alamut_server = "a-ht-na.interactive-biosoftware.com"
        String alamut_port = "80"
        String alamut_user_secret_name = "alamut-batch-ini-user"
        Int alamut_queue_limit = 4
        String alamut_queue_folder = "gs://biofx-task-queue/alamut"
        String alamut_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/alamut:20230630"
        Boolean alamut_save_working_files = false
        String alamut_anno_src_id = "228"
        String alamut_anno_min_age = "P6M"
		# FAST inputs
        Boolean has_haploid_sites = false
		String alamut_data_load_config_name = "Alamut"
    }

    # Pre-filter the collective/merged VCF to input to Alamut to remove
    #  variants that already have Alamut annotations
    call FASTUtils.FASTRemoveAlreadyAnnotatedFromVCFTask {
        input:
            input_vcf = merged_vcf,
            output_basename = batch_name + "_all.alamutprefilter",
            reference_build = reference_build,
            annotation_source_id = alamut_anno_src_id,
            annotation_min_age = alamut_anno_min_age,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    # Annotate with Alamut
    call AlamutBatch.AlamutBatchTask {
        input:
            input_vcf = FASTRemoveAlreadyAnnotatedFromVCFTask.output_vcf_gz,
            output_basename = batch_name + "_all.alamut",
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
    # Load Alamut annotations
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
    
    output {
        # Alamut vcf
        File alamut_vcf_gz = AlamutBatchTask.output_vcf_gz
        # Alamut data load ID
        String data_load_id = AlamutLoadTask.data_load_id
    }
}