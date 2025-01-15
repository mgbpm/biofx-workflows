version 1.0

import "../tasks/HelperTasks.wdl" as HelperTasks
import "../palantir/ScoringTasks.wdl" as ScoringTasks
import "../tasks/PRSStructs.wdl"

workflow PRSRawScoreWorkflow {
    input {
        File input_vcf
        File adjustment_model_manifest
    }

    AdjustmentModelData model_data = read_json(adjustment_model_manifest)

    # Clean up the query VCF
    call HelperTasks.RenameChromosomesInVcf as RenameVcf {
        input:
            vcf = input_vcf
    }

    call HelperTasks.GetBaseMemory as GetBaseMemoryFromVcf {
        input:
            vcf = RenameVcf.renamed
    }

    call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
        input:
            vcf = RenameVcf.renamed,
            mem = GetBaseMemoryFromVcf.gigabytes
    }

    # Check weights files and score VCF
    scatter (i in range(length(model_data.var_weights))) {
        File var_weights_file = model_data.var_weights[i]

        call ScoringTasks.CheckWeightsCoverSitesUsedInTraining {
            input:
                sites_used_in_training = model_data.training_variants,
                weight_set = object {linear_weights : var_weights_file}
        }

        call ScoringTasks.ScoreVcf {
            input:
                vcf = RenameVcf.renamed,
                basename = sub(basename(var_weights_file), ".var_weights.tsv", ""),
                weights = var_weights_file,
                base_mem = GetBaseMemoryFromVcf.gigabytes,
                sites = model_data.training_variants,
                chromosome_encoding = "MT"
        }
    }

    output {
        Array[File] prs_raw_scores = ScoreVcf.score
        Array[File] prs_raw_scores_log = ScoreVcf.log
        Array[File] sites_scored = ScoreVcf.sites_scored
    }
}
