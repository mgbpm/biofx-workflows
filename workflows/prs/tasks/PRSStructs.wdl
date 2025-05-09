version 1.0

# FOR INTERNAL USE ONLY
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
}
