version 1.0

import "PCATasks.wdl"
import "TrainAncestryAdjustmentModel.wdl"
import "Structs.wdl"
import "HelperTasks.wdl"

workflow MakeAdjustmentModelWorkflow {
  input {
    File   weights
    File   pca_variants
    File   reference_vcf
    String name
    File?  query_variants
  }

  call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInWeights {
    input:
        tsv        = weights
      , skipheader = true
  }

  call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInPcaVariants {
    input:
        tsv        = pca_variants
      , skipheader = false
  }

  call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInReferenceVcf {
    input:
      vcf = reference_vcf
  }

  call PCATasks.ArrayVcfToPlinkDataset as ReferenceBed {
    input:
        vcf             = RenameChromosomesInReferenceVcf.renamed
      , pruning_sites   = RenameChromosomesInPcaVariants.renamed
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
    linear_weights : RenameChromosomesInWeights.renamed
  }

  NamedWeightSet named_weight_set = object {
    condition_name : name,
    weight_set     : weight_set
  }

  call TrainAncestryAdjustmentModel.TrainAncestryAdjustmentModel as TrainModel {
    input:
        named_weight_set    = named_weight_set
      , population_pcs      = PerformPCA.pcs
      , population_vcf      = RenameChromosomesInReferenceVcf.renamed
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
