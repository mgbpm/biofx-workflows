version 1.0

import "PCATasks.wdl"
import "TrainAncestryAdjustmentModel.wdl"
import "Structs.wdl"

workflow MakeAdjustmentModelWorkflow {
  input {
    File   weights
    File   pca_variants
    File   reference_vcf
    String name
    File?  query_variants
  }

  call PCATasks.ArrayVcfToPlinkDataset as ReferenceBed {
    input:
        vcf             = reference_vcf
      , pruning_sites   = pca_variants
      , basename        = name
      , subset_to_sites = query_variants
  }

  call PCATasks.PerformPCA {
    input:
        bed      = ReferenceBed.bed
      , bim      = ReferenceBed.bim
      , fam      = ReferenceBed.fam
      , basename = name
  }

  WeightSet weight_set = object {
    linear_weights : weights
  }

  NamedWeightSet named_weight_set = object {
    condition_name : name,
    weight_set     : weight_set
  }

  call TrainAncestryAdjustmentModel.TrainAncestryAdjustmentModel as TrainModel {
    input:
        named_weight_set    = named_weight_set
      , population_pcs      = PerformPCA.pcs
      , population_vcf      = reference_vcf
      , population_basename = name
      , sites               = query_variants
  }

  output {
    File    parameters = TrainModel.fitted_params
    File    sites      = TrainModel.sites_used_in_scoring
    # .........................................................................
    File    pcs        = PerformPCA.pcs
    File    pcloadings = PerformPCA.pc_loadings
    File    pcmeansd   = PerformPCA.mean_sd

    # -------------------------------------------------------------------------

    Boolean converged  = TrainModel.fit_converged
    File    scores     = TrainModel.adjusted_population_scores
  }
}
