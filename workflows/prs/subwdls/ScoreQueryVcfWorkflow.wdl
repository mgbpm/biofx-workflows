version 1.0

import "../tasks/ScoringTasks.wdl"
import "../tasks/PCATasks.wdl"
import "../tasks/HelperTasks.wdl"
import "../tasks/Structs.wdl"

workflow ScoreQueryVcf {
  input {
    File    query_vcf
    File    adjustment_model_manifest
    String  name
    Boolean norename = false
  }

  AdjustmentModelData model_data = read_json(adjustment_model_manifest)

  if (! norename) {
    call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
      input:
          vcf = query_vcf
    }
  }

  File query_vcf_ = select_first([RenameChromosomesInQueryVcf.renamed,
                                  query_vcf])

  call HelperTasks.GetBaseMemory as GetBaseMemoryFromVcf {
    input:
        vcf = query_vcf_
  }

  Int base_memory = select_first([GetBaseMemoryFromVcf.gigabytes,
                                  model_data.base_memory])

  call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
    input:
        vcf = query_vcf_
      , mem = base_memory
  }

  # --------------------------------------------------------------------------

  WeightSet weight_set = object {
    linear_weights : model_data.scoring_inputs[0].variant_weights
  }

  call ScoringTasks.CheckWeightsCoverSitesUsedInTraining {
    input:
        sites_used_in_training = model_data.scoring_inputs[0].training_variants
      , weight_set             = weight_set
  }

  call ScoringTasks.ScoreVcf as ScoreQueryVcf {
    input:
        vcf                 = query_vcf_
      , weights             = model_data.scoring_inputs[0].variant_weights
      , sites               = model_data.scoring_inputs[0].training_variants
      , chromosome_encoding = "MT"
      , base_mem            = base_memory
      , basename            = name
  }

  call PCATasks.ArrayVcfToPlinkDataset as QueryBed {
    input:
        vcf           = query_vcf_
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
      , output_basename     = name
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
