version 1.0

# Forked from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

struct SelfExclusiveSites {
  File sites
  Int  maxAllowed
}

struct WeightSet {
  File                linear_weights
  File?               interaction_weights
  SelfExclusiveSites? interaction_self_exclusive_sites
}

struct ScoringInputs {
  File variant_weights
  File training_variants
}

struct AdjustmentModelData {
  String               condition_code
  File                 parameters
  File                 principal_components
  File                 loadings
  File                 meansd
  Array[File]          variant_weights
  File?                score_weights
  Array[ScoringInputs] scoring_inputs
  File                 pca_variants
  File                 original_pca_variants
  File                 query_file
  Int                  base_memory
  Float                population_mean
  Float                population_sd
}