version 1.0

import "ScoringTasks.wdl"
import "PCATasks.wdl"
import "HelperTasks.wdl"
import "Structs.wdl"

workflow ScoreVcf {
  input {
    File query_vcf
    File weights
    File pca_variants

    File model_parameters
    File training_variants
    File principal_components
    File loadings
    File meansd
  }

  call HelperTasks.GetBaseMemory {
    input:
      vcf = query_vcf
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

  call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
    input:
      vcf = query_vcf
  }

  WeightSet weight_set = object {
    linear_weights : RenameChromosomesInWeights.renamed
  }
  
  call ScoringTasks.CheckWeightsCoverSitesUsedInTraining {
    input:
        sites_used_in_training = training_variants
      , weight_set             = weight_set
  }

  call ScoringTasks.ScoreVcf as ScoreQueryVcf {
    input:
        vcf                 = RenameChromosomesInQueryVcf.renamed
      , weights             = RenameChromosomesInWeights.renamed
      , sites               = training_variants
      , chromosome_encoding = "MT"
      , base_mem            = GetBaseMemory.gigabytes
      , basename            = "query"
  }

  call PCATasks.ArrayVcfToPlinkDataset as QueryBed {
    input:
        vcf           = RenameChromosomesInQueryVcf.renamed
      , pruning_sites = RenameChromosomesInPcaVariants.renamed
      , mem           = GetBaseMemory.gigabytes
      , basename      = "temp"
  }

  call PCATasks.ProjectArray as ProjectQuery {
    input:
        pc_loadings = loadings
      , pc_meansd   = meansd
      , bed         = QueryBed.bed
      , bim         = QueryBed.bim
      , fam         = QueryBed.fam
      , mem         = GetBaseMemory.gigabytes
      , basename    = "query"
  }

  call ScoringTasks.AdjustScores {
    input:
        fitted_model_params = model_parameters
      , pcs                 = ProjectQuery.projections
      , scores              = ScoreQueryVcf.score
  }

  call ScoringTasks.MakePCAPlot {
    input:
        population_pcs = principal_components
      , target_pcs     = ProjectQuery.projections
  }

  output {
    File raw_scores      = ScoreQueryVcf.score
    File adjusted_scores = AdjustScores.adjusted_scores
    File pc_projection   = ProjectQuery.projections
    File pc_plot         = MakePCAPlot.pca_plot
  }
}
