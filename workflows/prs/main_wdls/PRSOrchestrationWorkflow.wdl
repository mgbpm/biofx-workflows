version 1.0

import "../../../steps/Utilities.wdl"
import "../../../steps/FileUtils.wdl"
import "../../lowpassimputation/Glimpse2Imputation.wdl"
import "../subwdls/PRSRawScoreWorkflow.wdl"
import "../subwdls/PRSMixScoreWorkflow.wdl"
import "../subwdls/PRSPCAWorkflow.wdl"
import "../tasks/PRSStructs.wdl"
import "../tasks/HelperTasks.wdl"

workflow PRSOrchestrationWorkflow {
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
        Boolean impute_reference_only_variants = false
        Boolean call_indels = false
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
        Boolean collect_glimpse_qc = true
        File? glimpse_monitoring_script
        File ref_fasta
        File ref_fai
        File ref_dict
        String glimpse_docker_image = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
        String glimpse_extract_docker_image = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"

        # PRS INPUTS
        String prs_test_code
        Array[File] model_manifests
        File conditions_config
        Boolean norename = false
        String ubuntu_docker_image = "ubuntu:latest"
        String prs_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250210"
    }

    # Fetch CRAM AND CRAI for sample
    call FileUtils.FetchFilesTask as FetchFiles {
        input:
            data_location = select_first([data_location]),
            file_types = ["cram", "crai"],
            recursive = false,
            file_match_keys = [subject_id, sample_id],
            docker_image = orchutils_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            disk_size = fetch_disk_size
    }

    # Run GLIMPSE to get imputed low-pass variants
    call Glimpse2Imputation.Glimpse2Imputation as RunGlimpse {
        input:
            reference_chunks = select_first([glimpse_reference_chunks]),
            crams = select_all(select_first([[FetchFiles.cram], []])),
            cram_indices = select_all(select_first([[FetchFiles.crai], []])),
            sample_ids = [sample_id],
            fasta = ref_fasta,
            fasta_index = ref_fai,
            ref_dict = ref_dict,
            output_basename = subject_id + "_" + sample_id + "_" + prs_test_code,
            impute_reference_only_variants = impute_reference_only_variants,
            call_indels = call_indels,
            n_burnin = n_burnin,
            n_main = n_main,
            effective_population_size = effective_population_size,
            collect_qc_metrics = collect_glimpse_qc,
            preemptible = 9,
            docker = glimpse_docker_image,
            docker_extract_num_sites_from_reference_chunk = glimpse_extract_docker_image,
            cpu_ligate = 4,
            mem_gb_ligate = 4,
            monitoring_script = glimpse_monitoring_script
    }

    scatter (i in range(length(model_manifests))) {
        String condition_name = sub(basename(model_manifests[i]), "_[[0-9]]+\\.json$", "")

        AdjustmentModelData model_data = read_json(model_manifests[i])

        # Get PRS raw scores for each condition
        call PRSRawScoreWorkflow.PRSRawScoreWorkflow as PRSRawScores {
            input:
                input_vcf = RunGlimpse.imputed_vcf,
                adjustment_model_manifest = model_manifests[i],
                norename = norename
        }

        # Get the PRS mix raw score for each condition
        call PRSMixScoreWorkflow.PRSMixScoreWorkflow as PRSMixScores {
            input:
                output_basename = condition_name,
                raw_scores = PRSRawScores.prs_raw_scores,
                score_weights = select_first([model_data.score_weights])
        }

        # Adjust the PRS mix raw score with PCA and model
        call PRSPCAWorkflow.PRSPCAWorkflow as PerformPCA {
            input:
                output_basename = condition_name,
                input_vcf = RunGlimpse.imputed_vcf,
                adjustment_model_manifest = model_manifests[i],
                prs_raw_scores = PRSMixScores.prs_mix_raw_score,
                norename = norename
        }
    }

    # Create summary of risk score, percentile, and condition info for reporting
    call SummarizeScores {
        input:
            condition_codes = select_first([select_all(condition_name)]),
            scores = select_first([select_all(PerformPCA.adjusted_scores)]),
            conditions_config = conditions_config,
            output_filename = subject_id + "_" + sample_id + "_" + prs_test_code + "_results.tsv",
            docker_image = prs_docker_image
    }

    output {
        # Glimpse outputs
        File glimpse_vcf = RunGlimpse.imputed_vcf
        File glimpse_vcf_index = RunGlimpse.imputed_vcf_index
        File? glimpse_qc_metrics = RunGlimpse.qc_metrics

        # PRS Outputs
        Array[Array[File]] prs_raw_scores = PRSRawScores.prs_raw_scores
        Array[File] prs_mix_raw_score = PRSMixScores.prs_mix_raw_score
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
        String output_filename
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
            --scores-list WORK/scores_files_list.txt \
            --condition-codes WORK/condition_codes_list.txt \
            --output-file "~{output_filename}.tsv" \
            --output-dir OUTPUT
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File risk_summary = "OUTPUT/~{output_filename}.tsv"
    }
}
