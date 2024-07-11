version 1.0

import "Structs.wdl"
import "ScoringPart.wdl"

struct RunSpec {
    File    imputed_array_vcf
    File    population_vcf
    File    pruning_sites_for_pca
    String  basename

    NamedWeightSet
            named_weight_set
}

workflow ScoringPartRunner {
  input {
    Array[RunSpec] specs
  }

  scatter (spec in specs) {
    call ScoringPart.ScoringImputedDataset {
      input:
        imputed_array_vcf     = spec.imputed_array_vcf    ,
        population_vcf        = spec.population_vcf       ,
        pruning_sites_for_pca = spec.pruning_sites_for_pca,
        basename              = spec.basename             ,
        named_weight_set      = spec.named_weight_set     ,

        redoPCA               = true,
        adjustScores          = true,
        population_pcs        = "PLACEHOLDER__REQUIRED_BY_BUGGY_CODE",
        population_meansd     = "PLACEHOLDER__REQUIRED_BY_BUGGY_CODE",
        population_loadings   = "PLACEHOLDER__REQUIRED_BY_BUGGY_CODE",
        population_basename   = "1000G",
    }
  }

  output {
    Array[File]  raw_scores        = ScoringImputedDataset.raw_scores
    Array[File?] pc_plots          = ScoringImputedDataset.pc_plot
    Array[File?] population_scores = ScoringImputedDataset.adjusted_population_scores
    Array[File?] imputed_scores    = ScoringImputedDataset.adjusted_array_scores
  }
}
