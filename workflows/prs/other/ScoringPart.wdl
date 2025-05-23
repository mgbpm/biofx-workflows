version 1.0

# Forked from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

import "../tasks/PCATasks.wdl" as PCATasks
import "../tasks/ScoringTasks.wdl" as ScoringTasks
import "TrainAncestryAdjustmentModel.wdl" as TrainAncestryAdjustmentModel
import "../tasks/Structs.wdl"

workflow ScoringImputedDataset {
	input {
	NamedWeightSet named_weight_set

	File imputed_array_vcf  # imputed VCF for scoring (and optionally PCA projection): make sure the variant IDs exactly match those in the weights file
	Int extract_ids_plink_mem = 8
	Int extract_ids_population_mem = 8
	Int scoring_mem = 16
	Int population_scoring_mem = scoring_mem * 4
	Int vcf_to_plink_mem = 8
	Int population_vcf_to_plink_mem = 8
	Int pca_memory = 8
	Int project_array_memory = 8


	String? population_basename # for naming the output of population scoring
	String basename # for naming the output of array scoring and the array projection files

	## these next 3 files are after performing PCA on your population dataset ( again, make sure all the variant IDs are the same )
	File? population_loadings
	File? population_meansd
	File? population_pcs
	File? pruning_sites_for_pca # and the sites used for PCA

	File? population_vcf # population VCF, output from PerformPopulationPCA.  The variant IDs must exactly match those in the weights file
	AncestryAdjustmentModelParams? fitted_model_params_and_sites # model parameters from fitting to reference population.  Either fitted_model_params or population_vcf must be passed as input

	String? columns_for_scoring # Plink expects the first 3 columns in your weights file to be variant ID, effect allele, effect weight
	# if this isn't true, then you should give it the correct column #s in that order
	# example: if you were to set columns_for_scoring = "11 12 13" would mean that the 11th column is the variant ID, the 12th column
	# is the effect allele, and the 13th column is the effect weight
	Boolean redoPCA = false
	Boolean adjustScores = true
	Boolean make_model_only = false
  }

  if (adjustScores) {
		#check for required optional inputs
		if (!defined(population_loadings)) {
			call ErrorWithMessage as ErrorPopulationLoadings {
				input:
					message = "Optional input population_loadings must be included if adjusting scores"
			}
		}

		if (!defined(population_meansd)) {
			call ErrorWithMessage as ErrorMeanSD {
				input:
					message = "Optional input population_meansd must be included if adjusting scores"
			}
		}

		if (!defined(population_pcs)) {
			call ErrorWithMessage as ErrorPopulationPCs {
				input:
					message = "Optional input population_pcs must be included if adjusting scores"
			}
		}

		if (!defined(pruning_sites_for_pca)) {
			call ErrorWithMessage as ErrorPruningSites {
				input:
					message = "Optional input pruning_sites_for_pca must be included if adjusting scores"
			}
	}

		if (!defined(fitted_model_params_and_sites) && !defined(population_vcf)) {
			call ErrorWithMessage as ErrorNeitherModelOrPopulation {
				input:
					message = "Either fitted_model_params or population_vcf must be passed as input if adjusting scores"
			}
		}

		if (defined(fitted_model_params_and_sites) && defined(population_vcf)) {
			call ErrorWithMessage as ErrorBothModelAndPopultion {
				input:
					message = "Only one of fitted_model_params and population_vcf should be passed as input if adjusting scores"
			}
		}
  }

	if (adjustScores && defined(population_vcf)) {
		call ScoringTasks.ExtractIDsPlink as ExtractIDsPopulation {
			input:
				vcf = select_first([population_vcf]),
				mem = extract_ids_population_mem
		}
	}

	if (defined(fitted_model_params_and_sites)) {
		AncestryAdjustmentModelParams local_params = select_first([fitted_model_params_and_sites])
		File sites_used_in_scoring_for_model = local_params.sites_used_in_scoring
		File fitted_params_for_model = local_params.fitted_model_params

		call ScoringTasks.CheckWeightsCoverSitesUsedInTraining {
		input:
			sites_used_in_training = sites_used_in_scoring_for_model,
			weight_set = named_weight_set.weight_set
		}
	}

	if (defined(ExtractIDsPopulation.ids) || defined(sites_used_in_scoring_for_model)) {
		File sites_to_use_in_scoring = select_first([ExtractIDsPopulation.ids, sites_used_in_scoring_for_model])
	}

        if (!make_model_only) {

            call ScoringTasks.DetermineChromosomeEncoding {
              input:
                  weights = named_weight_set.weight_set.linear_weights
            }

            call ScoringTasks.ScoreVcf as ScoreImputedArray {
              input:
                  vcf                 = imputed_array_vcf
                , basename            = basename
                , weights             = named_weight_set.weight_set.linear_weights
                , base_mem            = scoring_mem
                , extra_args          = columns_for_scoring
                , sites               = sites_to_use_in_scoring
                , chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding
            }

            if (defined(named_weight_set.weight_set.interaction_weights)) {
                call ScoringTasks.AddInteractionTermsToScore {
                  input:
                      vcf                  = imputed_array_vcf
                    , interaction_weights  = select_first([named_weight_set.weight_set.interaction_weights])
                    , scores               = ScoreImputedArray.score
                    , sites                = sites_to_use_in_scoring
                    , basename             = basename
                    , self_exclusive_sites = named_weight_set.weight_set.interaction_self_exclusive_sites
                }

                call ScoringTasks.CombineScoringSites {
                  input:
                      sites_used_linear_score      = ScoreImputedArray.sites_scored
                    , sites_used_interaction_score = AddInteractionTermsToScore.sites_used_in_interaction_score
                    , basename                     = basename
                }
            }
        }

	if (adjustScores) {
		call ScoringTasks.ExtractIDsPlink {
			input:
				vcf = imputed_array_vcf,
				mem = extract_ids_plink_mem
		}

		if (redoPCA && defined(population_vcf)) {
			call PCATasks.ArrayVcfToPlinkDataset as PopulationArrayVcfToPlinkDataset {
				input:
					vcf = select_first([population_vcf]),
					pruning_sites = select_first([pruning_sites_for_pca]),
					subset_to_sites = ExtractIDsPlink.ids,
					basename = "population",
					mem = population_vcf_to_plink_mem
			}

			call PCATasks.PerformPCA {
			  input:
				bim = PopulationArrayVcfToPlinkDataset.bim,
				bed = PopulationArrayVcfToPlinkDataset.bed,
				fam = PopulationArrayVcfToPlinkDataset.fam,
				mem = pca_memory,
				basename = basename
			}
		}

		call PCATasks.ArrayVcfToPlinkDataset {
			input:
			vcf = imputed_array_vcf,
			pruning_sites = select_first([pruning_sites_for_pca]),
			basename = basename,
			mem = vcf_to_plink_mem
		}

		if (defined(population_vcf)) {
			call CheckPopulationIdsValid {
				input:
					pop_vcf_ids = select_first([ExtractIDsPopulation.ids]),
					pop_pc_loadings = select_first([PerformPCA.pc_loadings, population_loadings]),
			}
		}

		call PCATasks.ProjectArray {
			input:
				pc_loadings = select_first([PerformPCA.pc_loadings, population_loadings]),
				pc_meansd = select_first([PerformPCA.mean_sd, population_meansd]),
				bed = ArrayVcfToPlinkDataset.bed,
				bim = ArrayVcfToPlinkDataset.bim,
				fam = ArrayVcfToPlinkDataset.fam,
				mem = project_array_memory,
				basename = basename
		}

		if (defined(population_vcf)) {
			call TrainAncestryAdjustmentModel.TrainAncestryAdjustmentModel {
				input:
					named_weight_set = named_weight_set,
					population_pcs = select_first([PerformPCA.pcs, population_pcs]),
					population_vcf = select_first([population_vcf]),
					population_basename = select_first([population_basename]),
					sites = ExtractIDsPlink.ids
			}
		}

                call ScoringTasks.MakePCAPlot {
                  input:
                      population_pcs = select_first([PerformPCA.pcs, population_pcs])
                    , target_pcs     = ProjectArray.projections
                }

                if (!make_model_only) {
                    call ScoringTasks.AdjustScores {
                      input:
                          fitted_model_params = select_first([TrainAncestryAdjustmentModel.fitted_params, fitted_params_for_model])
                        , pcs                 = ProjectArray.projections
                        , scores              = select_first([AddInteractionTermsToScore.scores_with_interactions, ScoreImputedArray.score])
                    }

                    call ScoringTasks.CompareScoredSitesToSitesUsedInTraining {
                      input:
                          sites_used_in_training = select_first([TrainAncestryAdjustmentModel.sites_used_in_scoring, sites_used_in_scoring_for_model])
                        , sites_used_in_scoring  = select_first([CombineScoringSites.combined_scoring_sites, ScoreImputedArray.sites_scored])
                        , weight_set             = named_weight_set.weight_set
                    }

                    if (CompareScoredSitesToSitesUsedInTraining.n_missing_sites > 0) {
                        #if there expected sites are missing, calculate potential effect on scores
                        call ScoringTasks.AddShiftToRawScores as ShiftScoresUpForMissingSites {
                          input:
                              raw_scores = select_first([AddInteractionTermsToScore.scores_with_interactions, ScoreImputedArray.score])
                            , shift      = CompareScoredSitesToSitesUsedInTraining.max_error_up
                            , basename   = "shifted_raw_scores_up_for_missing_sites"
                        }

                        call ScoringTasks.AddShiftToRawScores as ShiftScoresDownForMissingSites {
                          input:
                              raw_scores = select_first([AddInteractionTermsToScore.scores_with_interactions, ScoreImputedArray.score])
                            , shift      = CompareScoredSitesToSitesUsedInTraining.max_error_down
                            , basename   = "shifted_raw_scores_down_for_missing_sites"
                        }

                        call ScoringTasks.AdjustScores as AdjustScoresShiftedUpForMissingSites {
                          input:
                              fitted_model_params = select_first([TrainAncestryAdjustmentModel.fitted_params, fitted_params_for_model])
                            , pcs                 = ProjectArray.projections
                            , scores              = ShiftScoresUpForMissingSites.shifted_scores
                        }

                        call ScoringTasks.AdjustScores as AdjustScoresShiftedDownForMissingSites {
                          input:
                              fitted_model_params = select_first([TrainAncestryAdjustmentModel.fitted_params, fitted_params_for_model])
                            , pcs                 = ProjectArray.projections
                            , scores              = ShiftScoresDownForMissingSites.shifted_scores
                        }
                    }

                    call ScoringTasks.CombineMissingSitesAdjustedScores {
                      input:
                          adjusted_scores_shifted_up   = select_first([AdjustScoresShiftedUpForMissingSites.adjusted_scores, AdjustScores.adjusted_scores])
                        , adjusted_scores_shifted_down = select_first([AdjustScoresShiftedDownForMissingSites.adjusted_scores, AdjustScores.adjusted_scores])
                        , adjusted_scores              = AdjustScores.adjusted_scores
                        , n_missing_sites              = CompareScoredSitesToSitesUsedInTraining.n_missing_sites
                        , condition_name               = named_weight_set.condition_name
                    }
                }

		if (!select_first([CheckPopulationIdsValid.files_are_valid, true])) {
			call ErrorWithMessage {
				input:
				message = "Population VCF IDs are not a subset of the population PCA IDs; running with these inputs would give an incorrect result."
			}
		}
	}

  if (!make_model_only) {
    File raw_scores_ = select_first([AddInteractionTermsToScore.scores_with_interactions, ScoreImputedArray.score])
  }

  output {
    File?    pc_projection                 = ProjectArray.projections
    File?    raw_scores                    = raw_scores_
    File?    pc_plot                       = MakePCAPlot.pca_plot
    File?    adjusted_array_scores         = AdjustScores.adjusted_scores

    File?    model_parameters              = TrainAncestryAdjustmentModel.fitted_params
    File?    training_variants             = TrainAncestryAdjustmentModel.sites_used_in_scoring
    File?    pcs                           = PerformPCA.pcs
    File?    pcloadings                    = PerformPCA.pc_loadings
    File?    pcmeansd                      = PerformPCA.mean_sd

    Boolean? fit_converged                 = TrainAncestryAdjustmentModel.fit_converged
    File?    adjusted_population_scores    = TrainAncestryAdjustmentModel.adjusted_population_scores

    Int?     n_missing_sites_from_training = CompareScoredSitesToSitesUsedInTraining.n_missing_sites
    File?    missing_sites_shifted_scores  = CombineMissingSitesAdjustedScores.missing_sites_shifted_scores
  }
}


task CheckPopulationIdsValid{
	input {
		File pop_vcf_ids
		File pop_pc_loadings

	}
	command <<<
		# check if population VCF file contains a subset of population PC loading ids

		# 1. extract IDs, removing first rows of the pc files
		awk '{print $1}' ~{pop_pc_loadings} | tail -n +2 > pop__pc_ids.txt

		comm -23 <(sort pop_pc_ids.txt | uniq) <(sort ~{pop_vcf_ids} | uniq) > array_specific_ids.txt
		if [[ -s array_specific_ids.txt ]]
		then
		echo false
		else
		echo true
		fi

	>>>
	output {
		Boolean files_are_valid = read_boolean(stdout())
	}
	runtime {
		docker: "ubuntu:21.10"
	}
}

 #Print given message to stderr and return an error
 task ErrorWithMessage{
	 input {
		 String message
	 }
	command <<<
		>&2 echo "Error: ~{message}"
		exit 1
	>>>

	 runtime {
		 docker: "ubuntu:21.10"
	 }
 }
