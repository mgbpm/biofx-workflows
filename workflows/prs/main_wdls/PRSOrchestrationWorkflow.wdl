version 1.0

import "../../../steps/Utilities.wdl"
import "../../lowpassimputation/Glimpse2Imputation.wdl"
import "../subwdls/PRSRawScoreWorkflow.wdl"
import "../subwdls/PRSMixScoreWorkflow.wdl"
import "../subwdls/PRSPCAWorkflow.wdl"
import "../tasks/PRSStructs.wdl"
import "../tasks/HelperTasks.wdl"

workflow PRSOrchestrationWorkflow {
    input {
        # GLIMPSE2 INPUTS
        File? glimpse_reference_chunks
        File? input_cram
        File? input_crai
        String sample_id
        Boolean impute_reference_only_variants = false
        Boolean call_indels = false
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
        Boolean collect_glimpse_qc = true
        String glimpse_docker_image = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
        String glimpse_extract_docker_image = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
        File? glimpse_monitoring_script
        File? ref_fasta
        File? ref_fai
        File? ref_dict

        # PRS INPUTS
        String subject_id
        String prs_test_code
        Array[File] model_manifests
        File conditions_config
        String ubuntu_docker_image = "ubuntu:latest"
        String prs_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250210"
        Boolean norename = false
        
        # DEBUGGING INPUTS
        Boolean run_glimpse = true
        File? query_vcf
    }

    # Run input checks
    if (run_glimpse) {
        if (!defined(glimpse_reference_chunks) || !defined(input_cram) || !defined(input_crai)) {
            String GlimpseInputError = "Missing one or more GLIMPSE inputs: ref chunks, input cram, and/or input crai."
        }
        if (!defined(ref_fasta) || !defined(ref_fai) || !defined(ref_dict)) {
            String GlimpseReferenceInputError = "Missing one or more reference files for GLIMPSE: fasta, fai, and/or dict."
        }
    }

    if (!run_glimpse && !defined(query_vcf)) {
        String InputVcfError = "If GLIMPSE is not being run, please input a VCF for running PRS."
    }

    String input_check_result = select_first([GlimpseInputError, GlimpseReferenceInputError, InputVcfError, "No error"])

    if (input_check_result != "No error") {
        call Utilities.FailTask as FailInputCheck {
            input:
                error_message = input_check_result
        }
    }

    if (input_check_result == "No error") {
        if (run_glimpse) {
            # Run GLIMPSE to get imputed low-pass variants
            call Glimpse2Imputation.Glimpse2Imputation as RunGlimpse {
                input:
                    reference_chunks = select_first([glimpse_reference_chunks]),
                    crams = select_all(select_first([[input_cram], []])),
                    cram_indices = select_all(select_first([[input_crai], []])),
                    sample_ids = [sample_id],
                    fasta = select_first([ref_fasta]),
                    fasta_index = select_first([ref_fai]),
                    ref_dict = select_first([ref_dict]),
                    output_basename = subject_id + "_" + sample_id + "_" + prs_test_code,
                    impute_reference_only_variants = impute_reference_only_variants,
                    call_indels = call_indels,
                    n_burnin = select_first([n_burnin]),
                    n_main = select_first([n_main]),
                    effective_population_size = select_first([effective_population_size]),
                    collect_qc_metrics = collect_glimpse_qc,
                    preemptible = 9,
                    docker = glimpse_docker_image,
                    docker_extract_num_sites_from_reference_chunk = glimpse_extract_docker_image,
                    cpu_ligate = 4,
                    mem_gb_ligate = 4,
                    monitoring_script = select_first([glimpse_monitoring_script])
            }
        }

        scatter (model_manifest in model_manifests) {
            String condition_name = sub(basename(model_manifest), "_[[0-9]]+\\.json$", "")

            AdjustmentModelData model_data = read_json(model_manifest)

            # Get PRS raw scores for each condition
            call PRSRawScoreWorkflow.PRSRawScoreWorkflow as PRSRawScores {
                input:
                    query_vcf = select_first([RunGlimpse.imputed_vcf, query_vcf]),
                    adjustment_model_manifest = model_manifest,
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
                    input_vcf = select_first([RunGlimpse.imputed_vcf, query_vcf]),
                    adjustment_model_manifest = model_manifest,
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
    }

    output {
        # Glimpse outputs
        File? glimpse_vcf = RunGlimpse.imputed_vcf
        File? glimpse_vcf_index = RunGlimpse.imputed_vcf_index
        File? glimpse_qc_metrics = RunGlimpse.qc_metrics
        Array[File?]? glimpse_phase_monitoring = RunGlimpse.glimpse_phase_monitoring
        File? glimpse_ligate_monitoring = RunGlimpse.glimpse_ligate_monitoring

        # PRS Outputs
        Array[Array[File]]? prs_raw_scores = PRSRawScores.prs_raw_scores
        Array[File]? prs_mix_raw_score = PRSMixScores.prs_mix_raw_score
        Array[File?]? prs_adjusted_score = PerformPCA.adjusted_scores
        Array[File]? pc_projection = PerformPCA.pc_projection
        Array[File]? pc_plot = PerformPCA.pc_plot

        # Individual Outputs
        File? risk_summary = SummarizeScores.risk_summary
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
