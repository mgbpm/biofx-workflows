version 1.0

import "../lowpassimputation/Glimpse2Imputation.wdl"
import "PRSRawScoreWorkflow.wdl"
import "PRSMixScoreWorkflow.wdl"
import "PRSPCAWorkflow.wdl"
import "PRSAdjustmentWorkflow.wdl"
import "PRSCategorizationWorkflow.wdl"

workflow PRSOrchestrationWorkflow {
	input {
		# GLIMPSE2 INPUTS
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

		# PRS INPUTS
		Array[File] condition_zip_files
		File pruning_sites_for_pca
		String ubuntu_docker_image = "ubuntu:21.10"
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		String interaction_docker_image = "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e"
		String flash_pca_docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
		String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"

		# REFERENCE FILES
		File ref_fasta
		File ref_fai
		File ref_dict
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

	scatter (i in range(length(condition_zip_files))){
		call Unzip 
		# Get PRS raw scores for each condition
		call PRSRawScoreWorkflow.PRSRawScoreWorkflow as PRSRawScores {
			input:
				input_vcf = RunGlimpse.imputed_vcf,
				condition_file = condition_zip_files[i],
				plink_docker_image = plink_docker_image,
				interaction_docker_image = interaction_docker_image
		}
		# Get the PRS mix raw score for each condition
		call PRSMixScoreWorkflow.PRSMixScoreWorkflow as PRSMixScores {
			input:
				raw_scores = PRSRawScores.prs_raw_scores,
				condition_file = condition_zip_files[i],
				ubuntu_docker_image = ubuntu_docker_image
		}
		# Perform PCA with population model
		call PRSPCAWorkflow.PRSPCAWorkflow as PerformPCA {
			input:
				condition_file = condition_zip_files[i],
				input_vcf = RunGlimpse.imputed_vcf,
				pruning_sites_for_pca = pruning_sites_for_pca,
				chr_encoding = PRSRawScores.chromosome_encoding,
				plink_docker_image = plink_docker_image,
				flash_pca_docker_image = flash_pca_docker_image,
				tidyverse_docker_image = tidyverse_docker_image

		}
		# Adjust PRS mix score for each condition
		call PRSAdjustmentWorkflow.PRSAdjustmentWorkflow as AdjustPRSScores {
			input:
				condition_file = condition_zip_files[i],
				pca_projections = PerformPCA.pc_projection,
				prs_mix_raw_score = PRSMixScores.prs_mix_raw_score,
				tidyverse_docker_image = tidyverse_docker_image
		}
	}

	# Categorize each condition's score into bins; report percentile & bin
	call PRSCategorizationWorkflow.PRSCategorizationWorkflow as CategorizeScores {
		input:
			prs_scores = AdjustPRSScores.adjusted_scores,
			sample_ids = PRSMixScores.sample_ids[0],
			ubuntu_docker_image = ubuntu_docker_image
	}

	output {
		# Glimpse outputs
		File glimpse_vcf = RunGlimpse.imputed_vcf
        File glimpse_vcf_index = RunGlimpse.imputed_vcf_index
        File? glimpse_qc_metrics = RunGlimpse.qc_metrics
        Array[File?] glimpse_phase_monitoring = RunGlimpse.glimpse_phase_monitoring
        File? glimpse_ligate_monitoring = RunGlimpse.glimpse_ligate_monitoring

		# PRS Outputs
		File prs_mix_adjusted_score = AdjustPRSScores.adjusted_scores
		File pc_projection = PerformPCA.pc_projection
		File pc_plot = PerformPCA.pc_plot

		# Individual Outputs
		Array[File] individuals_risk_summaries = CategorizeScores.individual_risk_summaries
	}
}