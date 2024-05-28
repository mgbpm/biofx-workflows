version 1.0

import "ScoringTasks.wdl"
import "ScoreVcf.wdl"

workflow TrainModelWorkflow {
  input {
    File  vcf                   # population VCF
    File  linear_weights
    File  principal_components
    # File? scoring_variants      # aka imputed variants
  }

  call ScoreVcf.ScoreVcfWorkflow as Scoring {
    input:
        vcf      = vcf
      , weights  = linear_weights
      # , variants = scoring_variants  # aka imputed variants
  }

  call ScoringTasks.TrainAncestryModel as TrainModel {
    input:
        population_pcs    = principal_components
      , population_scores = Scoring.scores
      , output_basename   = "model"
  }

  output {
    Boolean converged  = TrainModel.fit_converged
    File    parameters = TrainModel.fitted_params
    File    variants   = Scoring.variants
    File    scores     = TrainModel.adjusted_population_scores
  }
}
