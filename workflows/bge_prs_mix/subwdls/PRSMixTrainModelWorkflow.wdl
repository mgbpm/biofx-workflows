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
		# Interaction inputs
		File? interaction_weights
		SelfExclusiveSites? self_exclusive_sites
		# Training model inputs
		File population_pcs
		# Docker images
		String python_docker_image = "python:3.9.10"
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		String ubuntu_docker_image = "ubuntu:21.10"
		String interaction_docker_image = "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e"
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

	String population_basename = sub(basename(population_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
	if (defined(interaction_weights)) {
		if (!defined(self_exclusive_sites)) {
			call Utilities.FailTask as SitesInputFail {
				input:
					error_message = "Need self exclusive sites object to add interaction terms to score."
			}
		}
		call ScoringTasks.AddInteractionTermsToScore as AddInteractionTermsToScorePopulation {
			input:
				vcf = population_vcf,
				interaction_weights = select_first([interaction_weights]),
				scores = GetMixScore.prs_mix_raw_score,
				sites = scoring_sites,
				basename = condition_name + "_" + population_basename,
				self_exclusive_sites = self_exclusive_sites,
				docker_image = interaction_docker_image
		}
		call ScoringTasks.CombineScoringSites {
			input:
				sites_used_linear_score = ScorePopulationVCF.prs_sites_scored[0],
				sites_used_interaction_score = AddInteractionTermsToScorePopulation.sites_used_in_interaction_score,
				basename = condition_name + "_" + population_basename,
				docker_image = ubuntu_docker_image
		}
	}

	call ScoringTasks.TrainAncestryModel {
		input:
			population_pcs = population_pcs,
			population_scores = select_first([AddInteractionTermsToScorePopulation.scores_with_interactions, ScorePopulationVCF.prs_raw_scores]),
			output_basename = condition_name + "_" + population_basename,
			docker_image = tidyverse_docker_image
	}

	output {
		File fitted_params = TrainAncestryModel.fitted_params
		File sites_used_in_scoring = select_first([CombineScoringSites.combined_scoring_sites, ScorePopulationVCF.prs_sites_scored[0]])
		File raw_population_scores = select_first([AddInteractionTermsToScorePopulation.scores_with_interactions, ScorePopulationVCF.prs_raw_scores])
		File adjusted_population_scores = TrainAncestryModel.adjusted_population_scores
		Boolean fit_converged = TrainAncestryModel.fit_converged
	}
}
