version 1.0

import "ScoringTasks.wdl"
import "PCATasks.wdl"
import "HelperTasks.wdl"
import "Structs.wdl"

workflow ScoreQueryVcf {
  input {
    File   query_vcf
    File   adjustment_model_manifest
    String name
  }

  AdjustmentModelData model_data = read_json(adjustment_model_manifest)

  call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
    input:
        vcf = query_vcf
  }

  call HelperTasks.GetBaseMemory as GetBaseMemoryFromVcf {
    input:
        vcf = RenameChromosomesInQueryVcf.renamed
  }

  Int  base_memory        = select_first([GetBaseMemoryFromVcf.gigabytes,
                                          model_data.base_memory])

  call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
    input:
        vcf = RenameChromosomesInQueryVcf.renamed
      , mem = base_memory
  }

  # --------------------------------------------------------------------------

  WeightSet weight_set = object {
    linear_weights : model_data.weights
  }

  call ScoringTasks.CheckWeightsCoverSitesUsedInTraining {
    input:
        sites_used_in_training = model_data.training_variants
      , weight_set             = weight_set
  }

  call ScoringTasks.ScoreVcf as ScoreQueryVcf {
    input:
        vcf                 = RenameChromosomesInQueryVcf.renamed
      , weights             = model_data.weights
      , sites               = model_data.training_variants
      , chromosome_encoding = "MT"
      , base_mem            = base_memory
      , basename            = name
  }

  call PCATasks.ArrayVcfToPlinkDataset as QueryBed {
    input:
        vcf           = RenameChromosomesInQueryVcf.renamed
      , pruning_sites = model_data.pca_variants
      , mem           = base_memory
      , basename      = "query"
  }

  call PCATasks.ProjectArray as ProjectQuery {
    input:
        pc_loadings = model_data.loadings
      , pc_meansd   = model_data.meansd
      , bed         = QueryBed.bed
      , bim         = QueryBed.bim
      , fam         = QueryBed.fam
      , mem         = base_memory
      , basename    = "query"
  }

  call ScoringTasks.AdjustScores {
    input:
        fitted_model_params = model_data.parameters
      , pcs                 = ProjectQuery.projections
      , scores              = ScoreQueryVcf.score
  }

  call ScoringTasks.MakePCAPlot {
    input:
        population_pcs = model_data.principal_components
      , target_pcs     = ProjectQuery.projections
  }

  output {
    File raw_scores      = ScoreQueryVcf.score
    File adjusted_scores = AdjustScores.adjusted_scores
    File pc_projection   = ProjectQuery.projections
    File pc_plot         = MakePCAPlot.pca_plot
    File kept_variants   = ExtractQueryVariants.ids
  }
}
