version 1.0

import "../lowpassimputation/Glimpse2Imputation.wdl"
import "PRSMixWorkflow.wdl"

workflow BGEPRSWorkflow {
	input {
		# Reference inputs
		File ref_fasta
		File ref_fai
		File ref_dict

		# Glimpse2 inputs
		File glimpse_reference_chunks
		Array[File] input_crams
        Array[File] input_crai
		Array[String] sample_ids
		Boolean impute_reference_only_variants = false
        Boolean call_indels = false
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
		Boolean collect_glimpse_qc = true
		String glimpse_docker_image = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
		String glimpse_extract_docker_image = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
		File? glimpse_monitoring_script

		# PRS inputs
		Array[File] var_weights
		File score_weight
		File population_vcf
		File pruning_sites_for_pca
	}

	# Run GLIMPSE to get imputed low-pass variants
	call Glimpse2Imputation.Glimpse2Imputation as RunGlimpse {
		input:
			reference_chunks = glimpse_reference_chunks,
			crams = input_crams,
			cram_indices = input_crai,
			sample_ids = sample_ids,
			fasta = ref_fasta,
			fasta_index = ref_fai,
			output_basename = "",
			ref_dict = ref_dict,
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

	# Run PRS Mix Workflow
	call PRSMixWorkflow as RunPRSMix {
		input:
			imputed_vcf = RunGlimpse.imputed_vcf,
			imputed_vcf_index = RunGlimpse.imputed_vcf_index
			var_weights = var_weights,
			score_weight = score_weight,
			population_vcf = population_vcf,
			pruning_sites_for_pca = pruning_sites_for_pca
	}

	# For each disease...
		# Bin individuals based on thresholds
		# Determine individual's percentile based on PRS score

	output {
		# Glimpse outputs
		File glimpse_vcf = RunGlimpse.imputed_vcf
        File glimpse_vcf_index = RunGlimpse.imputed_vcf_index
        File? glimpse_qc_metrics = RunGlimpse.qc_metrics
        Array[File?] glimpse_phase_monitoring = RunGlimpse.glimpse_phase_monitoring
        File? glimpse_ligate_monitoring = RunGlimpse.glimpse_ligate_monitoring

		# PRS Outputs
		File prs_mix_adjusted_score = RunPRSMix.adjusted_scores
		File pc_projection = RunPRSMix.pc_projection
		File pc_plot = RunPRSMix.pc_plot

		# Individual Outputs
		#Int/String individual_bin
		#Int individual_percentile
		#File individual_qc_file
	}
}