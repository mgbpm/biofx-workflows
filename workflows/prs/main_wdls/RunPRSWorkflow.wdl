version 1.0

import "../subwdls/MakeAdjustmentModelWorkflow.wdl"
import "../subwdls/ScoreQueryVcfWorkflow.wdl"
import "../tasks/Structs.wdl"

workflow RunPRS {

  input {
    File    weights
    File    pca_variants
    File    reference_vcf
    String  model_name
    File    query_vcf
    String? query_name
    Boolean norename      = false
  }

  call MakeAdjustmentModelWorkflow.MakeAdjustmentModel {
    input:
        weights       = [weights]
      , pca_variants  = pca_variants
      , reference_vcf = reference_vcf
      , query_file    = query_vcf
      , name          = model_name
      , norename      = norename
  }

  String resolved_query_name = select_first([query_name,
                                             basename(query_vcf, ".vcf.gz")])
  AdjustmentModelData model_data = read_json(MakeAdjustmentModel.adjustment_model_manifest)

  call ScoreQueryVcfWorkflow.ScoreQueryVcf {
    input:
        query_vcf                 = model_data.query_file
      , adjustment_model_manifest = MakeAdjustmentModel.adjustment_model_manifest
      , name                      = resolved_query_name
      , norename                  = true   # sic
  }

  output {
    File    adjustment_model_manifest = MakeAdjustmentModel.adjustment_model_manifest
    Boolean converged                 = MakeAdjustmentModel.converged
    File    raw_reference_scores      = MakeAdjustmentModel.raw_reference_scores[0]
    File    adjusted_reference_scores = MakeAdjustmentModel.adjusted_reference_scores

    # -------------------------------------------------------------------------

    File    raw_query_scores          = ScoreQueryVcf.raw_scores
    File    adjusted_query_scores     = ScoreQueryVcf.adjusted_scores
    File    pc_projection             = ScoreQueryVcf.pc_projection
    File    pc_plot                   = ScoreQueryVcf.pc_plot
  }
}

# -------------------------------------------------------------------------------
