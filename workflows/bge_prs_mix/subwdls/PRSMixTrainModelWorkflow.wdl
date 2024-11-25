version 1.0

# Adapted from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

import "../tasks/ScoringTasks.wdl"
import "../tasks/PRSStructs.wdl"
import "../../steps/Utilities.wdl"
import "PRSRawScoreWorkflow.wdl"
import "PRSMixScoreWorkflow.wdl"

workflow TrainPRSMixModel {
	input {
		# Scoring inputs
		String condition_name
		Array[File] var_weights
		File scoring_sites
		File population_vcf
		Int scoring_mem = 8
		File score_weights
		# Training model inputs
		File population_pcs
		# Docker images
		String python_docker_image = "python:3.9.10"
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		String ubuntu_docker_image = "ubuntu:21.10"
		String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
	}

	# Score population VCF with each var weights file
	call PRSRawScoreWorkflow.PRSRawScoreWorkflow as ScorePopulationVCF {
		input:
			var_weights = var_weights,
			scoring_sites = scoring_sites,
			input_vcf = population_vcf,
			scoring_mem = scoring_mem,
			python_docker_image = python_docker_image,
			plink_docker_image = plink_docker_image
	}

	# Calculate PRS mix score for population VCF
	call PRSMixScoreWorkflow.PRSMixScoreWorkflow as GetMixScore {
		input:
			condition_name = condition_name,
			raw_scores = ScorePopulationVCF.prs_raw_scores,
			score_weights = score_weights,
			ubuntu_docker_image = ubuntu_docker_image
	}

	call ScoringTasks.TrainAncestryModel {
		input:
			population_pcs = population_pcs,
			population_scores = GetMixScore.prs_mix_raw_score,
			output_basename = condition_name + "_MixModel",
			docker_image = tidyverse_docker_image
	}

	output {
		File fitted_params = TrainAncestryModel.fitted_params
		Array[File] raw_population_scores = ScorePopulationVCF.prs_raw_scores
		File mixed_population_scores = GetMixScore.prs_mix_raw_score
		File adjusted_population_scores = TrainAncestryModel.adjusted_population_scores
		Boolean fit_converged = TrainAncestryModel.fit_converged
	}
}
