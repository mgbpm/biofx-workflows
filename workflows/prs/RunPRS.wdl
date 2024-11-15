version 1.0

import "MakeAdjustmentModel.wdl"
import "ScoreVcfWorkflow.wdl"

workflow RunPRS {

  input {
    File    weights
    File    pca_variants
    File    reference_vcf
    String  model_name
    File    query_vcf
    String? query_name
  }

  call MakeAdjustmentModel.MakeAdjustmentModel {
    input:
        weights       = weights
      , pca_variants  = pca_variants
      , reference_vcf = reference_vcf
      , query_vcf     = query_vcf
      , name          = model_name
  }

  String resolved_query_name = select_first([query_name,
                                             basename(query_vcf, ".vcf.gz")])

  call ScoreVcfWorkflow.ScoreVcf {
    input:
        query_vcf                 = query_vcf
      , adjustment_model_manifest = MakeAdjustmentModel.adjustment_model_manifest
      , name                      = resolved_query_name
  }

  output {
    File    adjustment_model_manifest = MakeAdjustmentModel.adjustment_model_manifest
    Boolean converged                 = MakeAdjustmentModel.converged
    File    raw_reference_scores      = MakeAdjustmentModel.raw_reference_scores
    File    adjusted_reference_scores = MakeAdjustmentModel.adjusted_reference_scores

    # -------------------------------------------------------------------------

    File    raw_query_scores          = ScoreVcf.raw_scores
    File    adjusted_query_scores     = ScoreVcf.adjusted_scores
    File    pc_projection             = ScoreVcf.pc_projection
    File    pc_plot                   = ScoreVcf.pc_plot
  }
}

# -------------------------------------------------------------------------------
