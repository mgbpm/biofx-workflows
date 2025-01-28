version 1.0

# Adapted from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

import "../palantir/ScoringTasks.wdl" as ScoringTasks
import "../../../steps/Utilities.wdl"
import "PRSMixScoreWorkflow.wdl"

workflow PRSTrainMixModelWorkflow {
    input {
        # Scoring inputs
        String condition_name
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
      }

      if (nweights > 1) {
        # Calculate PRS mix score for population VCF
        call PRSMixScoreWorkflow.PRSMixScoreWorkflow as GetMixScore {
            input:
                condition_name = condition_name,
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
              output_basename = condition_name + "_MixModel"
      }

    }

    output {
        # Model outputs
        File fitted_params = select_first([TrainAncestryModel.fitted_params])
        Array[File] sites_used_in_scoring = select_first([ScorePopulationVCF.sites_scored])
        File adjusted_population_scores = select_first([TrainAncestryModel.adjusted_population_scores])
        Boolean fit_converged = select_first([TrainAncestryModel.fit_converged])
        # Scoring outputs
        Array[File] raw_population_scores = select_first([ScorePopulationVCF.score])
        File? mixed_population_scores = GetMixScore.prs_mix_raw_score
    }
}
