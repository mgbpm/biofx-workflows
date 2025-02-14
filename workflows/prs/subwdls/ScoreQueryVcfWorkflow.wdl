version 1.0

import "../tasks/PRSStructs.wdl"
import "../tasks/HelperTasks.wdl"
import "PRSPCAWorkflow.wdl"
import "PRSRawScoreWorkflow.wdl"

workflow ScoreQueryVcf {
  input {
    File    query_vcf
    File    adjustment_model_manifest
    String  name
    Boolean norename = false
  }

  if (! norename) {
    call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
      input:
          vcf = query_vcf
    }
  }

  File query_vcf_ = select_first([RenameChromosomesInQueryVcf.renamed,
                                  query_vcf])

  AdjustmentModelData model_data = read_json(adjustment_model_manifest)

  call PRSRawScoreWorkflow.PRSRawScoreWorkflow as RawScores {
    input:
        query_vcf                 = query_vcf_
      , adjustment_model_manifest = adjustment_model_manifest
      , norename                  = true
  }

  call PRSPCAWorkflow.PRSPCAWorkflow as PCA {
    input:
        output_basename           = name
      , input_vcf                 = query_vcf_
      , adjustment_model_manifest = adjustment_model_manifest
      , prs_raw_scores            = RawScores.prs_raw_scores[0]
      , norename                  = true
  }

  output {
    File raw_scores      = RawScores.prs_raw_scores[0]
    File adjusted_scores = select_first([PCA.adjusted_scores])
    File pc_projection   = PCA.pc_projection
    File pc_plot         = PCA.pc_plot
    File kept_variants   = RawScores.kept_variants
  }
}
