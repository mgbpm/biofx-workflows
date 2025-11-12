version 1.0

import "../tasks/HelperTasks.wdl" as HelperTasks
import "../tasks/ScoringTasks.wdl" as ScoringTasks
import "../tasks/Structs.wdl"

workflow RawScoreWorkflow {
    input {
        File input_vcf
        File adjustment_model_manifest
        Boolean norename = false
        File renaming_lookup = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
    }

    AdjustmentModelData model_data = read_json(adjustment_model_manifest)

    # Clean up the query VCF
    if (! norename) {
        call HelperTasks.RenameChromosomesInVcf as RenameVcf {
            input:
                vcf = input_vcf,
                rename = renaming_lookup
        }
    }

    File input_vcf_ = select_first([RenameVcf.renamed, input_vcf])

    call HelperTasks.GetBaseMemory as GetBaseMemoryFromVcf {
        input:
            vcf = input_vcf_
    }

    call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
        input:
            vcf = input_vcf_,
            mem = GetBaseMemoryFromVcf.gigabytes
    }

    # Check weights files and score VCF
    scatter (scoring_inputs in model_data.scoring_inputs) {

        File variant_weights = scoring_inputs.variant_weights

        call ScoringTasks.CheckWeightsCoverSitesUsedInTraining {
            input:
                sites_used_in_training = scoring_inputs.training_variants
              , weight_set             = object {
                                                  linear_weights : variant_weights
                                                }
        }

        call ScoringTasks.ScoreVcf {
            input:
                vcf                 = input_vcf_,
                basename            = basename(variant_weights),
                weights             = variant_weights,
                base_mem            = GetBaseMemoryFromVcf.gigabytes,
                sites               = scoring_inputs.training_variants,
                chromosome_encoding = "MT"
        }
    }

    output {
        Array[File] prs_raw_scores = ScoreVcf.score
        Array[File] prs_raw_scores_log = ScoreVcf.log
        Array[File] sites_scored = ScoreVcf.sites_scored
    }
}
