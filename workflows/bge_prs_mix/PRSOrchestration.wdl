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
		File? ref_fasta
		File? ref_fai
		File? ref_dict

		# PRS INPUTS
		Array[File] condition_zip_files
		File condition_yaml
		File? pruning_sites_for_pca
		String ubuntu_docker_image = "ubuntu:21.10"
		String python_docker_image = "python:3.11"
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		String interaction_docker_image = "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e"
		String flash_pca_docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
		String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"

		# DEBUGGING INPUTS
		Boolean run_glimpse = true
		Boolean run_scoring = true
		Boolean run_mix_scoring = true
		Boolean run_pca = true
		Boolean run_adjustment = true
		Boolean run_summary = true
		File? input_vcf
		File? pc_projections
		File? input_scores
	}

	# Run input checks
	if (run_glimpse) {
		if (!defined(glimpse_reference_chunks) || !defined(input_crams) || !defined(input_crai)) {
			call Utilities.FailTask as GlimpseInputFail {
				input:
					error_message = "Missing one or more GLIMPSE inputs: ref chunks, input crams, and/or input crai."
			}
		}
		if (!defined(ref_fasta) || !defined(ref_fai) || !defined(ref_dict)) {
			call Utilities.FailTask as GlimpseReferenceInputFail {
				input:
					error_message = "Missing one or more reference files for GLIMPSE: fasta, fai, and/or dict."
			}
		}
	}
	if (!run_glimpse && !defined(input_vcf)) {
		call Utilities.FailTask as InputVcfFail {
			input:
				error_message = "If GLIMPSE is not being run, please input a VCF for running PRS modules."
		}
	}


	if (run_glimpse) {
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
	}

	scatter (i in range(length(condition_zip_files))) {

		String condition_name = basename(condition_zip_files[i], ".tar")

		call UnzipConditionFile{
			input:
				condition_zip_file = condition_zip_files[i],
				basename = condition_name,
				docker_image = ubuntu_docker_image
		}

		if (run_scoring) {
			# Get PRS raw scores for each condition
			call PRSRawScoreWorkflow.PRSRawScoreWorkflow as PRSRawScores {
				input:
					condition_name = condition_name,
					var_weights = UnzipConditionFile.var_weights,
					scoring_sites = UnzipConditionFile.scoring_sites,
					input_vcf = select_first([RunGlimpse.imputed_vcf, input_vcf]),
					python_docker_image = python_docker_image
					plink_docker_image = plink_docker_image
			}
		}

		if (run_mix_scoring) {
			# Get the PRS mix raw score for each condition
			call PRSMixScoreWorkflow.PRSMixScoreWorkflow as PRSMixScores {
				input:
					condition_name = condition_name,
					raw_scores = PRSRawScores.prs_raw_scores,
					score_weights = UnzipConditionFile.score_weights,
					ubuntu_docker_image = ubuntu_docker_image
			}
		}

		if (run_pca) {
			# Perform PCA with population model
			call PRSPCAWorkflow.PRSPCAWorkflow as PerformPCA {
				input:
					condition_name = condition_name,
					var_weights = UnzipConditionFile.var_weights,
					input_vcf = select_first([RunGlimpse.imputed_vcf, input_vcf]),
					pc_loadings = UnzipConditionFile.pc_loadings,
					pc_meansd = UnzipConditionFile.pc_meansd,
					population_pcs = UnzipConditionFile.pcs,
					pruning_sites_for_pca = pruning_sites_for_pca,
					chr_encoding = select_first([PRSRawScores.chromosome_encoding]),
					python_docker_image = python_docker_image,
					plink_docker_image = plink_docker_image,
					flash_pca_docker_image = flash_pca_docker_image,
					tidyverse_docker_image = tidyverse_docker_image
			}
		}

		if (run_adjustment) {
			# Adjust PRS mix score for each condition
			call PRSAdjustmentWorkflow.PRSAdjustmentWorkflow as AdjustPRSScores {
				input:
					condition_name = condition_name,
					var_weights = UnzipConditionFile.var_weights,
					pca_projections = select_first([PerformPCA.pc_projection, pc_projections]),
					prs_mix_raw_score = select_first([PRSMixScores.prs_mix_raw_score, PRSRawScores.prs_raw_scores, input_scores]),
					fitted_model_params = UnzipConditionFile.fitted_model_params
					tidyverse_docker_image = tidyverse_docker_image
			}
		}
	}

	if (run_summary) {
		# Categorize each condition's score into bins; report percentile & bin
		call PRSSummaryWorkflow.PRSSummaryWorkflow as SummarizeScores {
			input:
				prs_scores = select_first([AdjustPRSScores.adjusted_scores, input_scores]),
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
		Array[File]? prs_raw_scores = PRSRawScores.prs_raw_scores
		File? prs_mix_raw_score = PRSMixScores.prs_mix_raw_score
		File? prs_adjusted_score = AdjustPRSScores.adjusted_scores
		File? pc_projection = PerformPCA.pc_projection
		File? pc_plot = PerformPCA.pc_plot

		# Individual Outputs
		Array[File]? individuals_risk_summaries = SummarizeScores.individual_risk_summaries
	}
}

task UnzipConditionFile {
	input {
		File condition_zip_file
		String basename
		String docker_image
		Int disk_size = ceil(size(condition_zip_file, "GB") * 2) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail
		tar -xf "~{condition_zip_file}"
	>>>

	runtime {
		docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + " GB"
        preemptible: preemptible
	}

	output {
		Array[File] var_weights = "~{basename}/~{basename}.harmonized_weights.txt"
		File fitted_model_params = "~{basename}/~{basename}.fitted_model_params.tsv"
		File pcs = "~{basename}/~{basename}.pc"
		File pc_loadings = "~{basename}/~{basename}.pc.loadings"
		File pc_meansd = "~{basename}/~{basename}.pc.meansd"
		File score_weights = "~{basename}/~{basename}.score_weights.txt"
		File scoring_sites = "~{basename}/~{basename}.sscore.vars"
	}
}