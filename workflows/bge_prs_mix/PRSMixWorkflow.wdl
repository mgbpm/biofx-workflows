version 1.0

import "PRSTasks.wdl"

workflow PRSMixWorkflow {
	input {
		# Scoring inputs
		File input_vcf
		File input_vcf_idx
		Array[File] var_weights_files
		# PRS Mix inputs
		File score_weights_file
		# Condition-specific adjustment inputs
		File population_loadings
		File population_meansd
		File population_pcs
		File ancestry_adjustment_model
		# Other adjustment inputs
		File pruning_sites_for_pca
		File scoring_sites
		# Docker images
		String ubuntu_docker_image = "ubuntu:21.10"
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		String python_docker_image = "python:3.9.10"
		String flash_pca_docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
		String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
	}

	String condition_name = sub(basename(score_weights_file), ".score_weights.txt", "")

	# Calculate raw PRS score per variant weights file
	scatter (i in range(length(var_weights_files))) {
		call PRSTasks.DetermineChromosomeEncoding as ChrEncoding {
			input:
				weights = var_weights_files[i],
				output_basename = condition_name
				docker_image = python_docker_image
		}
		call PRSTasks.ScoreVCF as RawScore {
			input:
				input_vcf = input_vcf,
				var_weights = var_weights_files[i],
				chromosome_encoding = ChrEncoding.chromosome_encoding,
				sites = scoring_sites,
				output_basename = condition_name,
				docker_image = plink_docker_image
		}
	}

	# Calculate raw PRS mix score
	call PRSTasks.CalculateMixScore as MixRawScore {
		input:
			raw_scores = RawScore.score,
			score_weights = score_weights_file,
			output_basename = condition_name
	}

	# Calculate PCA for individual
	call PRSTasks.ExtractIDs as {
		input:
			input_vcf = input_vcf,
			chromosome_encoding = ChrEncoding.chromosome_encoding
	}
	call PRSTasks.ArrayVCFToPlinkDataset as VCFToPlinkDataset {
		input:
		input_vcf = input_vcf,
		pruning_sites = pruning_sites_for_pca,
		chromosome_encoding = ChrEncoding.chromosome_encoding,
		output_basename = condition_name,
		docker_image = plink_docker_image
	}
	call PRSTasks.ProjectArray as ProjectPC {
		input:
			pc_loadings = population_loadings
			pc_meansd = population_meansd,
			bed = VCFToPlinkDataset.bed,
			bim = VCFToPlinkDataset.bim,
			fam = VCFToPlinkDataset.fam,
			output_basename = condition_name + "_pca",
			docker_image = flash_pca_docker_image
	}
	call PRSTasks.MakePCAPlot as PCAPlot {
		input:
			population_pcs = population_pcs,
			target_pcs = ProjectPC.projections,
			output_basename = condition_name,
			docker_image = tidyverse_docker_image
	}

	# Adjust score with model and PCA
	call PRSTasks.AdjustScores as AdjustedScores {
		input:
			model_parameters = ancestry_adjustment_model,
			pcs = ProjectPC.projections,
			scores = MixRawScore.prs_mix_raw_score,
			output_basename = condition_name,
			docker_image = tidyverse_docker_image
	}
	
	output {
		# Sample IDs list
		File sample_ids_list = MixRawScore.sample_ids_list
		# PRS output
		Array[File] prs_raw_scores = RawScore.score
		Array[File] prs_raw_scores_log = RawScore.log
    	Array[File] prs_raw_sites_scored = RawScore.sites_scored
		File prs_mix_raw_score = MixRawScore.prs_mix_raw_score
		File adjusted_scores = AdjustedScores.adjusted_scores
		# PCA output
		File pc_projection = ProjectPC.projections
		File pc_plot = PCAPlot.pca_plot
	}
}