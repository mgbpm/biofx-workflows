version 1.0

# Adapted from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

import "../tasks/ScoringTasks.wdl" as ScoringTasks
import "../../../steps/Utilities.wdl"
import "MixScoreWorkflow.wdl"

workflow TrainMixModelWorkflow {
    input {
        # Scoring inputs
        String condition_code
        Array[File]+ var_weights
        File scoring_sites
        File reference_vcf
        File? score_weights
        Int scoring_mem
        # Training model inputs
        File population_pcs
    }

    Int nweights = length(var_weights)

    if (nweights > 1 && !defined(score_weights)) {
      call Utilities.FailTask as MissingScoreWeights {
        input:
          error_message = "missing score_weights file (it is required when multiple variant weight are specified)"
      }
    }

    if (nweights == 1 || defined(score_weights)) {

      # Score reference VCF
      scatter (weights_file in var_weights) {
        call ScoringTasks.ScoreVcf as ScorePopulationVCF {
          input:
            vcf = reference_vcf,
            basename = basename(weights_file, ".tsv"),
            weights = weights_file,
            base_mem = scoring_mem,
            sites = scoring_sites,
            chromosome_encoding = "MT"
        }

        # The prepending of a empty string to the values of the Object
        # below effectively coerces their types from File to String.
        Object scoring_inputs_ = object {
            variant_weights   : "" + weights_file
          , training_variants : "" + ScorePopulationVCF.sites_scored
        }
      }

      if (nweights > 1) {
        # Calculate PRS mix score for population VCF
        call MixScoreWorkflow.MixScoreWorkflow as GetMixScore {
            input:
                output_basename = condition_code + "_" + basename(reference_vcf),
                raw_scores = ScorePopulationVCF.score,
                score_weights = select_first([score_weights]),
        }
      }

      File population_scores = select_first([GetMixScore.prs_mix_raw_score,
                                             ScorePopulationVCF.score[0]])

      # Train the ancestry adjustment model
      call ScoringTasks.TrainAncestryModel {
          input:
              population_pcs = population_pcs,
              population_scores = population_scores,
              output_basename = condition_code + "_MixModel"
      }

    }

    output {
        # Model outputs
        File           fitted_params              = select_first([TrainAncestryModel.fitted_params])
        Array[File]    sites_used_in_scoring      = select_first([ScorePopulationVCF.sites_scored])
        File           adjusted_population_scores = select_first([TrainAncestryModel.adjusted_population_scores])
        Boolean        fit_converged              = select_first([TrainAncestryModel.fit_converged])
        Array[Object]+ scoring_inputs             = select_first([scoring_inputs_])
        Array[File]+   variant_weights            = var_weights
        # Scoring outputs
        Array[File]    raw_population_scores      = select_first([ScorePopulationVCF.score])
        File?          mixed_population_scores    = GetMixScore.prs_mix_raw_score
    }
}
