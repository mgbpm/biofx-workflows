version 1.0

# WDL Script to run just the QC Eval subcomponent of the BGWGS Pipeline

import "../../steps/VCFUtils.wdl"
import "../../steps/QCEval.wdl"
import "../../steps/Utilities.wdl"
import "../../steps/FileUtils.wdl"

workflow QCEvalSubComponent {
    input {
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
        # bcftools docker image
        String bcftools_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        # subject, sample id and data location
        String subject_id
        String sample_id
        String sample_data_location

        Array[String] fetch_vcf_filter_keys = [subject_id, sample_id]
        Array[FileMatcher]? fetch_vcf_file_matchers
        Boolean fetch_files_verbose = false
        Int fetch_disk_size = 75

       # reference genome files
        String reference_build = "GRCh38"
        File ref_dict
        File ref_fasta
        File ref_fasta_index

        # vcf filter inputs
        File target_roi_bed = "gs://lmm-reference-data/roi/targetROI_hg38_2023_08_24_withCHR.bed"
        # qceval inputs
        String qceval_project_type = "BGE_DRAGEN_TP_BINNING"
        String qceval_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20260521"
        File? thresholds = "gs://lmm-reference-data/annotation/pmeval/thresholds_20250912.tsv"
        File? difficult_to_map_regions = "gs://lmm-reference-data/annotation/pmeval/difficult_to_map_regions_20250912.tgz"
    }

    # Fetch the VCF files from the storage location
    call FileUtils.FetchFilesTask as FetchVcf {
        input:
            data_location = sample_data_location,
            file_types = if defined(fetch_vcf_file_matchers) then [] else ["vcf"],
            recursive = false,
            file_match_keys = if defined(fetch_vcf_file_matchers) then [] else fetch_vcf_filter_keys,
            file_matchers = fetch_vcf_file_matchers,
            verbose = fetch_files_verbose,
            docker_image = orchutils_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name
    }

    File sample_vcf = select_first([FetchVcf.vcf])
    
    # Filter called VCF to target region
    call VCFUtils.FilterVCFWithBEDTask {
        input:
            input_vcf = sample_vcf,
            input_bed = target_roi_bed,
            output_basename = subject_id + "_" + sample_id + ".target",
            docker_image = "~{bcftools_docker_image}"
    }


    # The assignments at the end of the conditional block below ensure
    # that `ref_fasta` and `ref_fasta_index` are passed to
    # `QCEval.QCEvalTask` only when needed, avoiding unnecessary
    # localization--especially of `ref_fasta`, which is typically
    # large.  (Because these inputs have type `File`, they are
    # mandatory and cannot be omitted.  If passed directly to the
    # task, they would be localized unconditionally.  Assigning them
    # via a conditional block adds them to the task's inputs only when
    # needed.  In contrast, the other inputs specific to the
    # "BGE_DRAGEN_TP_BINNING" case--`thresholds` and
    # `difficult_to_map_regions`--have type `File?`, so they can
    # simply be left out when not needed; localization is not an issue
    # in this case.)

    if (qceval_project_type == "BGE_DRAGEN_TP_BINNING") {
      if (!(defined(thresholds) && defined(difficult_to_map_regions))) {
        call Utilities.FailTask as MissingQcevalBgeDragenTpBinningInputs {
            input:
                error_message = "Some inputs required for QCEval with qceval_project_type = \\\"~{qceval_project_type}\\\" are missing"
        }
      }
      File maybe_reference_fasta = ref_fasta
      File maybe_reference_fasta_fai = ref_fasta_index
    }

    # Annotate target VCF with QCEval INFO field
    call QCEval.QCEvalTask {
        input:
            input_vcf = FilterVCFWithBEDTask.output_vcf_gz,
            project_type = qceval_project_type,
            output_basename = subject_id + "_" + sample_id + ".qceval",
            reference_fasta = maybe_reference_fasta,
            reference_fasta_fai = maybe_reference_fasta_fai,
            thresholds_tsv = thresholds,
            regions_tgz = difficult_to_map_regions,
            docker_image = qceval_docker_image
    }

    output {
        # filtered VCF
        File target_vcf_gz = FilterVCFWithBEDTask.output_vcf_gz
        # annotated VCFs
        File qceval_vcf_gz = QCEvalTask.output_vcf_gz
    }
}