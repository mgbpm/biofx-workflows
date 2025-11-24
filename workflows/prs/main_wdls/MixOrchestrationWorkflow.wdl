version 1.0

import "../../../steps/Utilities.wdl"
import "../../../steps/FileUtils.wdl"
import "../../lowpassimputation/Glimpse2Imputation.wdl"
import "../subwdls/RunPRSmixWorkflow.wdl"
import "../tasks/Structs.wdl"
import "../tasks/HelperTasks.wdl"

workflow MixOrchestrationWorkflow {
    input {
        # FETCH FILE INPUTS
        String  data_location
        String  sample_id
        String  subject_id
        String  orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20250203"
        String  gcp_project_id         = "mgb-lmm-gcp-infrast-1651079146"
        String  workspace_name
        Int     fetch_disk_size        = 75

        # GLIMPSE2 INPUTS
        File    glimpse_reference_chunks
        File    ref_fasta
        File    ref_fai
        File    ref_dict
        File    gnomadAF_ref_vcf
        
        # PRS INPUTS
        String             prs_test_code
        Array[Array[File]] model_manifests
        File               conditions_config
        File               percentiles
        File?              score_weights
        Boolean            mix_before_adjustment = true
        Boolean            norename              = false
        File               renaming_lookup       = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
        String             ubuntu_docker_image   = "ubuntu:latest"
        String             prs_docker_image      = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515"
    }

    call FileUtils.FetchFilesTask as FetchFiles {
        input:
            data_location = data_location,
            file_types = ["cram", "crai"],
            recursive = false,
            file_match_keys = [subject_id, sample_id],
            docker_image = orchutils_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            disk_size = fetch_disk_size
    }

    call Glimpse2Imputation.Glimpse2Imputation as RunGlimpse {
        input:
            reference_chunks = glimpse_reference_chunks,
            crams = select_all([FetchFiles.cram]),
            cram_indices = select_all([FetchFiles.crai]),
            sample_ids = [sample_id],
            fasta = ref_fasta,
            fasta_index = ref_fai,
            ref_dict = ref_dict,
            output_basename = subject_id + "_" + sample_id + "_" + prs_test_code,
            af_cutoff = ">=0.0001",
            gnomadAF_ref_vcf = gnomadAF_ref_vcf,
    }

    scatter (condition_models in model_manifests) {
        AdjustmentModelData model_data = read_json(condition_models[0])
        String condition_code = model_data.condition_code

        call RunPRSmixWorkflow.RunPrsMixWorkflow as GetMixScores {
            input:
                query_vcf = RunGlimpse.imputed_afFiltered_vcf,
                score_weights = score_weights,
                adjustment_model_manifests = condition_models,
                condition_code = condition_code,
                mix_before_adjustment = mix_before_adjustment,
                norename = norename,
                ubuntu_docker_image = ubuntu_docker_image
        }
    }

    call SummarizeScores {
        input:
            condition_codes = condition_code,
            scores = GetMixScores.final_score,
            conditions_config = conditions_config,
            percentiles = percentiles,
            basename = subject_id + "_" + sample_id + "_" + prs_test_code + "_results",
            docker_image = prs_docker_image
    }

    output {
        File        glimpse_vcf        = RunGlimpse.imputed_afFiltered_vcf
        File?       glimpse_qc_metrics = RunGlimpse.qc_metrics
        Array[File] prs_score          = GetMixScores.final_score
        File        risk_summary       = SummarizeScores.risk_summary
    }
}

task SummarizeScores {
    input {
        Array[String] condition_codes
        Array[File] scores
        File conditions_config
        File percentiles
        String basename
        String docker_image
        Int disk_size = ceil(size(scores, "GB") + size(conditions_config, "GB")) + 10
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        mkdir -p WORK
        mkdir -p OUTPUT

        # Make a list of scores files for python to access
        for s in '~{sep="' '" scores}'; do
            printf -- "${s}\n" >> WORK/scores_files_list.txt
        done

        # Make a list of condition_names
        for c in '~{sep="' '" condition_codes}'; do
            printf -- "${c}\n" >> WORK/condition_codes_list.txt
        done

        summarize_mix_scores.py \
            --condition-config "~{conditions_config}" \
            --percentiles "~{percentiles}" \
            --scores-list WORK/scores_files_list.txt \
            --condition-codes WORK/condition_codes_list.txt \
            --output-file "~{basename}.tsv" \
            --output-dir OUTPUT
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File risk_summary = "OUTPUT/~{basename}.tsv"
    }
}
