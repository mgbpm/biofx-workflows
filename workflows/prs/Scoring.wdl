version 1.0

import "PCATasks.wdl"
import "ScoringTasks.wdl"
import "TrainAncestryAdjustmentModel.wdl"
import "Structs.wdl"

struct PcaResults {
  File loadings
  File meansd
  File pc
}

workflow ScoringImputedDataset {
  input {

    NamedWeightSet
             named_weight_set

    File     imputed_array_vcf               # imputed VCF for scoring
                                             # (and optionally PCA
                                             # projection): make sure
                                             # the variant IDs exactly
                                             # match those in the
                                             # weights file

    String   basename                        # for naming the output
                                             # of array scoring and
                                             # the array projection
                                             # files

    Int      scoring_mem            = 16
    Int      vcf_to_plink_mem       = 8

    Boolean  redoPCA                = false
    Boolean  adjustScores           = true

    PcaResults?
             population_pca

    File?    pca_variants           # the sites used for PCA

    # -------------------------------------------------------------------------

    File?    population_vcf                  # population VCF, output
                                             # from
                                             # PerformPopulationPCA.
                                             # The variant IDs must
                                             # exactly match those in
                                             # the weights file.

    AncestryAdjustmentModelParams?
             model   # model parameters from
                                             # fitting to reference
                                             # population.  Either
                                             # fitted_model_params or
                                             # population_vcf must be
                                             # passed as input

    # -------------------------------------------------------------------------

    String?  columns_for_scoring             # plink expects the first
                                             # 3 columns in your
                                             # weights file to be
                                             # variant ID, effect
                                             # allele, effect weight
                                             # if this isn't true,
                                             # then you should give it
                                             # the correct column
                                             # numbers in that order
                                             # example: if you were to
                                             # set columns_for_scoring
                                             # = "11 12 13" would mean
                                             # that the 11th column is
                                             # the variant ID, the
                                             # 12th column is the
                                             # effect allele, and the
                                             # 13th column is the
                                             # effect weight
  }

  # ---------------------------------------------------------------------------

  if (adjustScores) {

    Boolean pcaInputsOk = (   defined(population_pca)
                           || (defined(population_vcf) && redoPCA))

    if (!pcaInputsOk) {
      call ErrorWithMessage as MissingPcaInputs {
        input:
          message =   "If adjustScores is true, then either "
                    + "population_pca must be specified, or else "
                    + "population_vcf must be specified and redoPCA "
                    + "must be true."
      }
    }

    if (pcaInputsOk) {
      Boolean pruningInputOk = defined(pca_variants)

      if (!pruningInputOk) {
        call ErrorWithMessage as MissingPcaSites {
          input:
            message =   "Input pca_variants must be specified if "
                      + "adjustScores is true."
        }
      }
    }

    if (select_first([pruningInputOk, false])) {
      Boolean fittedInputOk = (defined(model) != defined(population_vcf))

      if (!fittedInputOk) {
        call ErrorWithMessage as NotModelXorPopulation {
          input:
            message =   "Exactly one of fitted_model_params or population_vcf "
                      + "must be specified if adjustScores is true."
        }
      }

    }

  } #$#

  Boolean inputsOk = select_first([fittedInputOk, false]) || ! adjustScores

  if (inputsOk) {

    if (defined(population_pca)) {
      PcaResults just_pca = select_first([population_pca])
      File loadings = just_pca.loadings
      File meansd   = just_pca.meansd
      File pcs      = just_pca.pc
    }

    if (adjustScores && defined(population_vcf)) {
      call ScoringTasks.ExtractIDsPlink as PopulationVariants {
        input:
          vcf = select_first([population_vcf])
      }
    }

    if (defined(model)) {
      AncestryAdjustmentModelParams
          local_params                     = select_first([model])
      File sites_used_in_scoring_for_model = local_params.sites_used_in_scoring
      File fitted_params_for_model         = local_params.fitted_model_params

      call ScoringTasks.CheckWeightsCoverSitesUsedInTraining {
        input:
          sites_used_in_training = sites_used_in_scoring_for_model,
          weight_set             = named_weight_set.weight_set
      }
    } #$#

    if (   defined(PopulationVariants.ids)
        || defined(sites_used_in_scoring_for_model)) {
      File sites_to_use_in_scoring =
          select_first([PopulationVariants.ids,
                        sites_used_in_scoring_for_model])
    } #$#

    call ScoringTasks.DetermineChromosomeEncoding {
      input:
        weights = named_weight_set.weight_set.linear_weights
    } #$#

    call ScoringTasks.ScoreVcf as ScoreImputedArray {
      input:
        vcf                 = imputed_array_vcf,
        basename            = basename,
        weights             = named_weight_set.weight_set.linear_weights,
        base_mem            = scoring_mem,
        extra_args          = columns_for_scoring,
        sites               = sites_to_use_in_scoring,
        chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding
    } #$#

    if (defined(named_weight_set.weight_set.interaction_weights)) {
      call ScoringTasks.AddInteractionTermsToScore {
        input:
          vcf                  = imputed_array_vcf,
          interaction_weights  = select_first([named_weight_set
                                               .weight_set
                                               .interaction_weights]),
          scores               = ScoreImputedArray.score,
          sites                = sites_to_use_in_scoring,
          basename             = basename,
          self_exclusive_sites = named_weight_set.weight_set
                                 .interaction_self_exclusive_sites
      }

      call ScoringTasks.CombineScoringSites {
        input:
          sites_used_linear_score      = ScoreImputedArray.sites_scored,
          sites_used_interaction_score = AddInteractionTermsToScore
                                         .sites_used_in_interaction_score,
          basename                     = basename
      }
    } #$#

    if (adjustScores) {

      call ScoringTasks.ExtractIDsPlink as ImputedVariants {
        input:
          vcf = imputed_array_vcf
      }

      if (redoPCA) {

        call PCATasks.ArrayVcfToPlinkDataset as PopulationBed {
          input:
            vcf             = select_first([population_vcf]),
            pruning_sites   = select_first([pca_variants]),
            basename        = "population",
            mem             = vcf_to_plink_mem,
            subset_to_sites = ImputedVariants.ids
        }

        call PCATasks.PerformPCA {
          input:
            bim      = PopulationBed.bim,
            bed      = PopulationBed.bed,
            fam      = PopulationBed.fam,
            basename = basename
        }

      }

      call PCATasks.ArrayVcfToPlinkDataset as ImputedBed {
        input:
          vcf           = imputed_array_vcf,
          pruning_sites = select_first([pca_variants]),
          basename      = basename,
          mem           = vcf_to_plink_mem
      }

      if (defined(population_vcf)) {
        call CheckPopulationIdsValid {
          input:
            pop_vcf_ids     = select_first([PopulationVariants.ids]),
            pop_pc_loadings = select_first([PerformPCA.pc_loadings,
                                            loadings]),
        }
      }

      call PCATasks.ProjectArray {
        input:
          pc_loadings = select_first([PerformPCA.pc_loadings, loadings]),
          pc_meansd   = select_first([PerformPCA.mean_sd,     meansd]),
          bed         = ImputedBed.bed,
          bim         = ImputedBed.bim,
          fam         = ImputedBed.fam,
          basename    = basename
      }

      if (defined(population_vcf)) {
        call TrainAncestryAdjustmentModel.TrainAncestryAdjustmentModel {
          input:
            named_weight_set    = named_weight_set,
            population_pcs      = select_first([PerformPCA.pcs, pcs]),
            population_vcf      = select_first([population_vcf]),
            population_basename = "population_basename",
            sites               = ImputedVariants.ids
        }
      }

      call ScoringTasks.AdjustScores {
        input:
          fitted_model_params = select_first([TrainAncestryAdjustmentModel
                                              .fitted_params,
                                              fitted_params_for_model]),
          pcs                 = ProjectArray.projections,
          scores              = select_first([AddInteractionTermsToScore
                                              .scores_with_interactions,
                                              ScoreImputedArray.score])
      }

      call ScoringTasks.MakePCAPlot {
        input:
          population_pcs = select_first([PerformPCA.pcs, pcs]),
          target_pcs     = ProjectArray.projections
      }

      call ScoringTasks.CompareScoredSitesToSitesUsedInTraining {
        input:
          sites_used_in_training = select_first([
                                       TrainAncestryAdjustmentModel
                                       .sites_used_in_scoring,
                                       sites_used_in_scoring_for_model]),
          sites_used_in_scoring = select_first([
                                       CombineScoringSites
                                       .combined_scoring_sites,
                                       ScoreImputedArray.sites_scored]),
          weight_set = named_weight_set.weight_set
      }


      if (CompareScoredSitesToSitesUsedInTraining.n_missing_sites > 0) {

        # if there expected sites are missing, calculate potential
        # effect on scores

        call ScoringTasks.AddShiftToRawScores as ShiftScoresUpForMissingSites {
          input:
            raw_scores = select_first([AddInteractionTermsToScore
                                       .scores_with_interactions,
                                       ScoreImputedArray.score]),
            shift      = CompareScoredSitesToSitesUsedInTraining.max_error_up,
            basename   = "shifted_raw_scores_up_for_missing_sites"
        }

        call ScoringTasks.AddShiftToRawScores
             as ShiftScoresDownForMissingSites {
          input:
            raw_scores = select_first([AddInteractionTermsToScore
                                       .scores_with_interactions,
                                       ScoreImputedArray.score]),
            shift      = CompareScoredSitesToSitesUsedInTraining.max_error_down,
            basename   = "shifted_raw_scores_down_for_missing_sites"
        }

        call ScoringTasks.AdjustScores
             as AdjustScoresShiftedUpForMissingSites {
          input:
            fitted_model_params = select_first([TrainAncestryAdjustmentModel
                                                .fitted_params,
                                                fitted_params_for_model]),
            pcs                 = ProjectArray.projections,
            scores              = ShiftScoresUpForMissingSites.shifted_scores
        }

        call ScoringTasks.AdjustScores
             as AdjustScoresShiftedDownForMissingSites {
          input:
            fitted_model_params = select_first([TrainAncestryAdjustmentModel
                                                .fitted_params,
                                                fitted_params_for_model]),
            pcs                 = ProjectArray.projections,
            scores              = ShiftScoresDownForMissingSites.shifted_scores
        }
      }


      call ScoringTasks.CombineMissingSitesAdjustedScores {
        input:

          condition_name               = named_weight_set.condition_name,
          adjusted_scores              = AdjustScores.adjusted_scores,
          n_missing_sites              =
              CompareScoredSitesToSitesUsedInTraining.n_missing_sites,
          adjusted_scores_shifted_up   =
              select_first([
                  AdjustScoresShiftedUpForMissingSites
                  .adjusted_scores,
                  AdjustScores.adjusted_scores]),
          adjusted_scores_shifted_down =
              select_first([
                  AdjustScoresShiftedDownForMissingSites
                  .adjusted_scores,
                  AdjustScores.adjusted_scores])
      }

      if (!select_first([CheckPopulationIdsValid.files_are_valid, true])) {
        call ErrorWithMessage as InvalidPcaIds {
          input:
            message =   "Population VCF IDs are not a superset of the "
                      + "population PCA IDs; running with these inputs would "
                      + "lead to incorrect results."
        }
      }
    }
  } #$#

  # ---------------------------------------------------------------------------
  
  output {
  
    File     raw_scores                    = select_first([
                                                 AddInteractionTermsToScore
                                                 .scores_with_interactions,
                                                 ScoreImputedArray.score])
  
    File?    pc_projection                 = ProjectArray.projections
  
    File?    pc_plot                       = MakePCAPlot.pca_plot
  
    File?    adjusted_population_scores    = TrainAncestryAdjustmentModel
                                             .adjusted_population_scores
  
    File?    adjusted_array_scores         = AdjustScores.adjusted_scores
  
    Boolean? fit_converged                 = TrainAncestryAdjustmentModel
                                             .fit_converged
  
    Int?     n_missing_sites_from_training = CompareScoredSitesToSitesUsedInTraining
                                             .n_missing_sites
  
    File?    missing_sites_shifted_scores  = CombineMissingSitesAdjustedScores
                                             .missing_sites_shifted_scores

    File     s_variants                    = select_first([
                                                 CombineScoringSites.combined_scoring_sites,
                                                 ScoreImputedArray.sites_scored
                                             ])

    File     t_variants                    = select_first([
                                                 TrainAncestryAdjustmentModel.sites_used_in_scoring,
                                                 sites_used_in_scoring_for_model
                                             ])
  }
}

task CheckPopulationIdsValid {
  input {
    File pop_vcf_ids
    File pop_pc_loadings
  }

  String stdout = 'STDOUT'

  command <<<

  #     comm -23 <(sort pop_pc_ids.txt | uniq) <(sort ~{pop_vcf_ids} | uniq) > array_specific_ids.txt
  #
  # The command above will create a NON-EMPTY file array_specific_ids.txt IFF
  # the file pop_pc_ids.txt includes IDs that are NOT included in the file
  # ~{pop_vcf_ids}, i.e. IFF the set of IDs mentioned in the file
  # pop_pc_ids.txt IS NOT A SUBSET of the set of IDs mentioned in the file
  # ~{pop_vcf_ids}.

  # check whether the set of variant IDs in the population VCF file is a
  # (proper) subset of the set of population PCA loadings ids

  cut                    \
      --fields=1         \
      ~{pop_pc_loadings} \
    | tail --lines=+2    \
    | sort --unique      \
    > pca_ids

  comm -2 -3 pca_ids <(sort ~{pop_vcf_ids}) > pca_ids_not_in_vcf

  if [[ -s pca_ids_not_in_vcf ]]
  then
      RESULT=false
  else
      RESULT=true
  fi

  printf -- '%s' "${RESULT}" > '~{stdout}'
  >>>

  output {
    Boolean files_are_valid = read_boolean(stdout)
  }

  runtime {
    docker: "ubuntu:21.10"
  }
}

task ErrorWithMessage{
  input {
    String message
  }

  command <<<
  printf -- 'Error: ~{message}\n' >&2
  exit 1
  >>>

  runtime {
    docker: "ubuntu:21.10"
  }
}
