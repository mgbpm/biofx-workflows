version 1.0

import "../../../steps/Utilities.wdl"
import "../../../steps/FileUtils.wdl"
import "../../lowpassimputation/Glimpse2Imputation.wdl"
import "../subwdls/RawScoreWorkflow.wdl"
import "../subwdls/MixScoreWorkflow.wdl"
import "../subwdls/AdjustScoreWorkflow.wdl"
import "../tasks/Structs.wdl"
import "../tasks/HelperTasks.wdl"

workflow MixOrchestrationWorkflow {
    input {
        # FETCH FILE INPUTS
        String data_location
        String sample_id
        String subject_id
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20250203"
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        Int fetch_disk_size = 75

        # GLIMPSE2 INPUTS
        File glimpse_reference_chunks
        File ref_fasta
        File ref_fai
        File ref_dict
        String glimpse_af_cutoff = ">=0.0001"
        File gnomadAF_ref_vcf
        Boolean impute_reference_only_variants = false
        Boolean call_indels = false
        Boolean keep_monomorphic_ref_sites = false
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
        Boolean collect_glimpse_qc = true
        Int glimpse_preemptible = 1
        File? glimpse_monitoring_script
        String glimpse_docker_image = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
        String glimpse_extract_docker_image = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"

        # PRS INPUTS
        String prs_test_code
        Array[File] model_manifests
        File conditions_config
        File percentiles
        Boolean norename = false
        File renaming_lookup = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
        String ubuntu_docker_image = "ubuntu:latest"
        String prs_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515"
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
            af_cutoff = glimpse_af_cutoff,
            gnomadAF_ref_vcf = gnomadAF_ref_vcf,
            impute_reference_only_variants = impute_reference_only_variants,
            call_indels = call_indels,
            keep_monomorphic_ref_sites = keep_monomorphic_ref_sites,
            n_burnin = n_burnin,
            n_main = n_main,
            effective_population_size = effective_population_size,
            collect_qc_metrics = collect_glimpse_qc,
            preemptible = glimpse_preemptible,
            docker = glimpse_docker_image,
            docker_extract_num_sites_from_reference_chunk = glimpse_extract_docker_image,
            monitoring_script = glimpse_monitoring_script
    }

    scatter (i in range(length(model_manifests))) {
        AdjustmentModelData model_data = read_json(model_manifests[i])

        String condition_code = model_data.condition_code

        call RawScoreWorkflow.RawScoreWorkflow as RawScores {
            input:
                input_vcf = RunGlimpse.imputed_afFiltered_vcf,
                adjustment_model_manifest = model_manifests[i],
                norename = norename,
                renaming_lookup = renaming_lookup
        }

        call MixScoreWorkflow.MixScoreWorkflow as MixScores {
            input:
                output_basename = condition_code,
                raw_scores = RawScores.prs_raw_scores,
                score_weights = select_first([model_data.score_weights])
        }

        call AdjustScoreWorkflow.AdjustScoreWorkflow as PerformPCA {
            input:
                output_basename = condition_code,
                input_vcf = RunGlimpse.imputed_afFiltered_vcf,
                adjustment_model_manifest = model_manifests[i],
                prs_raw_scores = MixScores.prs_mix_raw_score,
                norename = norename,
                renaming_lookup = renaming_lookup
        }
    }

    call SummarizeScores {
        input:
            condition_codes = condition_code,
            scores = select_first([select_all(PerformPCA.adjusted_scores)]),
            conditions_config = conditions_config,
            percentiles = percentiles,
            basename = subject_id + "_" + sample_id + "_" + prs_test_code + "_results",
            docker_image = prs_docker_image
    }

    output {
        # Glimpse outputs
        File glimpse_vcf = RunGlimpse.imputed_afFiltered_vcf
        File glimpse_vcf_index = RunGlimpse.imputed_afFiltered_vcf_index
        File? glimpse_qc_metrics = RunGlimpse.qc_metrics

        # PRS Outputs
        Array[Array[File]] prs_raw_scores = RawScores.prs_raw_scores
        Array[File] prs_mix_raw_score = MixScores.prs_mix_raw_score
        Array[File?] prs_adjusted_score = PerformPCA.adjusted_scores

        # Individual Outputs
        File risk_summary = SummarizeScores.risk_summary
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
