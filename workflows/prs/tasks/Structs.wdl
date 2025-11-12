version 1.0

# Forked from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

struct SelfExclusiveSites {
	File sites # must have columns id, chrom, pos
	Int maxAllowed
}

struct WeightSet {
	File linear_weights # standard prs weights file
	File? interaction_weights # interaction weights file, must have columns id_1, id_2, chrom_1, chrom_2, pos_1, pos_2, allele_1, allele_2, weight (order not important)
	SelfExclusiveSites? interaction_self_exclusive_sites # The interaction term will only be added in no more than selfExclusiveSites.maxAllowed of the
																												# effect alleles listed in SelfExclusizeSites.sites is observed
}

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