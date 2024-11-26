version 1.0

import "../../steps/Utilities.wdl"
import "../lowpassimputation/Glimpse2Imputation.wdl"
import "PRSRawScoreWorkflow.wdl"
import "PRSMixScoreWorkflow.wdl"
import "PRSPCAWorkflow.wdl"
import "PRSAdjustmentWorkflow.wdl"
import "PRSSummaryWorkflow.wdl"

workflow PRSOrchestrationWorkflow {
    input {
        # GLIMPSE2 INPUTS
        File? glimpse_reference_chunks
        Array[File]? input_crams
        Array[File]? input_crai
        Array[String]? sample_ids
        String? vcf_basename
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
        Array[File] condition_jsons
        File condition_yaml
        File pruning_sites_for_pca
        Int prs_scoring_mem = 8
        String ubuntu_docker_image = "ubuntu:21.10"
        String python_docker_image = "python:3.11"
        String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        String interaction_docker_image = "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e"
        String flash_pca_docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
        String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"

        # DEBUGGING INPUTS
        Boolean run_glimpse = true
        File? query_vcf
    }

    # Run input checks
    if (run_glimpse) {
        if (!defined(glimpse_reference_chunks) || !defined(input_crams) || !defined(input_crai)) {
            String GlimpseInputError = "Missing one or more GLIMPSE inputs: ref chunks, input crams, and/or input crai."
        }
        if (!defined(ref_fasta) || !defined(ref_fai) || !defined(ref_dict)) {
            String GlimpseReferenceInputError = "Missing one or more reference files for GLIMPSE: fasta, fai, and/or dict."
        }
    }

    if (!run_glimpse && !defined(query_vcf)) {
        String InputVcfError = "If GLIMPSE is not being run, please input a VCF for running PRS modules."
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
                    crams = select_first([input_crams]),
                    cram_indices = select_first([input_crai]),
                    sample_ids = select_first([sample_ids]),
                    fasta = select_first([ref_fasta]),
                    fasta_index = select_first([ref_fai]),
                    ref_dict = select_first([ref_dict]),
                    output_basename = select_first([vcf_basename, "PRS_Glimpse"]),
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

        scatter (i in range(length(condition_jsons))) {
            String condition_name = sub(basename(condition_jsons[i]), "_[[0-9]]+\\.json$", "")
            AdjustmentModelData model_data = read_json(adjustment_model_manifest)

            # Get PRS raw scores for each condition
            call PRSRawScoreWorkflow.PRSRawScoreWorkflow as PRSRawScores {
                input:
                    var_weights = model_data.var_weights,
                    scoring_sites = model_data.training_variants,
                    query_vcf = select_first([RunGlimpse.imputed_vcf, query_vcf]),
                    scoring_mem = prs_scoring_mem,
                    plink_docker_image = plink_docker_image
            }

            # Get the PRS mix raw score for each condition
            call PRSMixScoreWorkflow.PRSMixScoreWorkflow as PRSMixScores {
                input:
                    condition_name = condition_name,
                    raw_scores = PRSRawScores.prs_raw_scores,
                    score_weights = model_data.score_weights,
                    ubuntu_docker_image = ubuntu_docker_image
            }

            # Perform PCA with population model
            call PRSPCAWorkflow.PRSPCAWorkflow as PerformPCA {
                input:
                    condition_name = condition_name,
                    query_vcf = select_first([RunGlimpse.imputed_vcf, query_vcf]),
                    pc_loadings = model_data.loadings,
                    pc_meansd = model_data.meansd,
                    population_pcs = model_data.principal_components,
                    pruning_sites_for_pca = pruning_sites_for_pca,
                    weights_chr_encoding = PRSRawScores.chromosome_encoding[0],
                    plink_docker_image = plink_docker_image,
                    flash_pca_docker_image = flash_pca_docker_image,
                    tidyverse_docker_image = tidyverse_docker_image
            }

            # Adjust PRS mix score for each condition
            call PRSAdjustmentWorkflow.PRSAdjustmentWorkflow as AdjustPRSScores {
                input:
                    condition_name = condition_name,
                    weights_chr_encoding = PRSRawScores.chromosome_encoding[0],
                    pca_projections = PerformPCA.pc_projection,
                    prs_raw_scores = PRSMixScores.prs_mix_raw_score,
                    fitted_model_params = model_data.parameters,
                    tidyverse_docker_image = tidyverse_docker_image
            }
        }

        # Categorize each condition's score into bins; report percentile & bin
        call PRSSummaryWorkflow.PRSSummaryWorkflow as SummarizeScores {
            input:
                prs_scores = AdjustPRSScores.adjusted_scores,
                condition_yaml = condition_yaml,
                python_docker_image = python_docker_image
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
        Array[File]? prs_adjusted_score = AdjustPRSScores.adjusted_scores
        Array[File]? pc_projection = PerformPCA.pc_projection
        Array[File]? pc_plot = PerformPCA.pc_plot

        # Individual Outputs
        Array[File]? individuals_risk_summaries = SummarizeScores.individual_risk_summaries
    }
}
