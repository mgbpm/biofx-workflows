version 1.0

import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/ScoringTasks.wdl" as ScoringTasks

workflow PRSMixWorkflow {
	input {
		# Scoring inputs
		File imputed_vcf
		File imputed_vcf_idx
		Array[File] var_weights
		# PRS Mix inputs
		File score_weight
		# Adjustment inputs
		File population_vcf
		File pruning_sites_for_pca
	}

	# Scatter over the files within the zipped file and run scoring on each
	scatter (i in range(length(var_weights))) {
		String output_basename = sub(basename(var_weights[i]), "txt", "")
		call ScoringTasks.DetermineChromosomeEncoding as DetermineChrEncoding {
			input:
				weights = var_weights[i]
		}
		call ScoringTasks.ScoreVcf as ScoreImputedVCF {
			input:
				vcf = imputed_vcf,
				basename = output_basename,
				weights = var_weights[i],
				base_mem = 16,
				chromosome_encoding = DetermineChrEncoding.chromosome_encoding,
		}
	}

	# Calculate raw PRS Mix score
	call CalculateMixScore {
		input:
			raw_scores = ScoreImputedVCF.score,
			score_weight = score_weight
	}

	# Adjust raw PRS Mix score
	call ScoringTasks.ExtractIDsPlink {
		input:
			vcf = imputed_vcf,
			chromosome_encoding = DetermineChrEncoding.chromosome_encoding,
			mem = 8
	}
	call PCATasks.ArrayVcfToPlinkDataset as PopulationArrayVcfToPlinkDataset {
		input:
			vcf = population_vcf,
			pruning_sites = pruning_sites_for_pca,
			subset_to_sites = ExtractIDsPlink.ids,
			basename = "population",
			use_ref_alt_for_ids = use_ref_alt_for_ids,
			chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding
	}
	call PCATasks.PerformPCA {
	  input:
		bim = PopulationArrayVcfToPlinkDataset.bim,
		bed = PopulationArrayVcfToPlinkDataset.bed,
		fam = PopulationArrayVcfToPlinkDataset.fam,
		basename = basename
	}
	call PCATasks.ArrayVcfToPlinkDataset {
		input:
		vcf = imputed_array_vcf,
		pruning_sites = pruning_sites_for_pca,
		basename = basename,
		mem = 8,
		use_ref_alt_for_ids = use_ref_alt_for_ids,
		chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding
	}
	call CheckPopulationIdsValid {
		input:
			pop_vcf_ids = ExtractIDsPopulation.ids,
			pop_pc_loadings = PerformPCA.pc_loadings
	}
	call PCATasks.ProjectArray {
		input:
			pc_loadings = PerformPCA.pc_loadings
			pc_meansd = PerformPCA.mean_sd,
			bed = ArrayVcfToPlinkDataset.bed,
			bim = ArrayVcfToPlinkDataset.bim,
			fam = ArrayVcfToPlinkDataset.fam,
			basename = "pca"
	}
	call ScoringTasks.AdjustScores {
		input:
			fitted_model_params = fitted_params_for_model,
			pcs = ProjectArray.projections,
			scores = ScoreImputedVCF.score
	}
	call ScoringTasks.MakePCAPlot {
		input:
			population_pcs = PerformPCA.pcs,
			target_pcs = ProjectArray.projections
	}
	call ScoringTasks.CompareScoredSitesToSitesUsedInTraining {
		input:
			sites_used_in_training = sites_used_in_scoring_for_model,
			sites_used_in_scoring = ScoreImputedVCF.sites_scored,
			weight_set = named_weight_set.weight_set
	}
	if (CompareScoredSitesToSitesUsedInTraining.n_missing_sites > 0) {
		#if there expected sites are missing, calculate potential effect on scores
		call ScoringTasks.AddShiftToRawScores as ShiftScoresUpForMissingSites {
			input:
				raw_scores = select_first([AddInteractionTermsToScore.scores_with_interactions, prs_scores]),
				shift = CompareScoredSitesToSitesUsedInTraining.max_error_up,
				basename = "shifted_raw_scores_up_for_missing_sites"
		}

		call ScoringTasks.AddShiftToRawScores as ShiftScoresDownForMissingSites {
			input:
				raw_scores = select_first([AddInteractionTermsToScore.scores_with_interactions, prs_scores]),
				shift = CompareScoredSitesToSitesUsedInTraining.max_error_down,
				basename = "shifted_raw_scores_down_for_missing_sites"
		}

		call ScoringTasks.AdjustScores as AdjustScoresShiftedUpForMissingSites {
			input:
				fitted_model_params = select_first([TrainAncestryAdjustmentModel.fitted_params, fitted_params_for_model]),
				pcs = ProjectArray.projections,
				scores = ShiftScoresUpForMissingSites.shifted_scores
		}

		call ScoringTasks.AdjustScores as AdjustScoresShiftedDownForMissingSites {
			input:
				fitted_model_params = select_first([TrainAncestryAdjustmentModel.fitted_params, fitted_params_for_model]),
				pcs = ProjectArray.projections,
				scores = ShiftScoresDownForMissingSites.shifted_scores
		}
	}

	call ScoringTasks.CombineMissingSitesAdjustedScores {
		input:
			adjusted_scores_shifted_up = select_first([AdjustScoresShiftedUpForMissingSites.adjusted_scores, AdjustScores.adjusted_scores]),
			adjusted_scores_shifted_down = select_first([AdjustScoresShiftedDownForMissingSites.adjusted_scores, AdjustScores.adjusted_scores]),
			adjusted_scores = AdjustScores.adjusted_scores,
			n_missing_sites = CompareScoredSitesToSitesUsedInTraining.n_missing_sites,
			condition_name = named_weight_set.condition_name
	}

	output {
		# Scores
		Array[File] prs_raw_scores = ScoreImputedVCF.score
		Array[File] prs_raw_scores_log = ScoreImputedVCF.log
    	Array[File] prs_raw_sites_scored = ScoreImputedVCF.sites_scored
		File prs_mix_raw_score = CalculateMixScore.prs_mix_raw_score
		File adjusted_scores = AdjustScores.adjusted_scores
		# PCA
		File pc_projection = ProjectArray.projections
		File pc_plot = MakePCAPlot.pca_plot
	}
}

task CalculateMixScore {
	input {
		Array[File] raw_scores
		File score_weight
		String docker_image
		Int disk_size = ceil(size(score_weights_file, "GB")) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<


	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		File prs_mix_raw_score = "temp"
	}
}