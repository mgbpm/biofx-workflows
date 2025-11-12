version 1.0

import "../tasks/HelperTasks.wdl" as HelperTasks
import "../tasks/ScoringTasks.wdl" as ScoringTasks
import "../tasks/Structs.wdl"

workflow RawScoreWorkflow {
    input {
        File    input_vcf
        File    adjustment_model_manifest
        Boolean norename        = false
        File    renaming_lookup = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
    }

    AdjustmentModelData model_data = read_json(adjustment_model_manifest)

    if (!norename) {
        call HelperTasks.RenameChromosomesInVcf as RenameVcf {
            input:
                vcf = input_vcf,
                rename = renaming_lookup
        }
    }

    File input_vcf_ = select_first([RenameVcf.renamed, input_vcf])

    call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
        input:
            vcf = input_vcf_,
            mem = model_data.base_memory
    }

    scatter (scoring_inputs in model_data.scoring_inputs) {
        File variant_weights = scoring_inputs.variant_weights
        File training_variants = scoring_inputs.training_variants

        call ScoringTasks.CheckWeightsCoverSitesUsedInTraining {
            input:
                sites_used_in_training = training_variants,
                weight_set = object {linear_weights : variant_weights}
        }

        call ScoringTasks.ScoreVcf {
            input:
                vcf = input_vcf_,
                basename = basename(variant_weights),
                weights = variant_weights,
                base_mem = model_data.base_memory,
                sites = training_variants,
                chromosome_encoding = "MT"
        }
    }

    output {
        File        renamed_vcf  = input_vcf_
        Array[File] raw_scores   = ScoreVcf.score
        Array[File] sites_scored = ScoreVcf.sites_scored
    }
}
