version 1.0

import "../tasks/HelperTasks.wdl"
import "../tasks/ScoringTasks.wdl"

workflow PRSRawScoreWorkflow {
    input {
        # PRS inputs
        File input_vcf
        File adjustment_model_manifest
        # Docker images
        String python_docker_image = "python:3.9.10"
        String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
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
                chromosome_encoding = "MT",
                sites = model_data.training_variants,
                weights = var_weights_file,
                basename = sub(basename(var_weights_file), ".var_weights.tsv", ""),
                base_mem = GetBaseMemoryFromVcf.gigabytes,
                docker_image = plink_docker_image
        }
    }

    output {
        # Chromosme encoding
        Array[String] chromosome_encoding = ChrEncoding.chromosome_encoding
        # PRS outputs
        Array[File] prs_raw_scores = ScoreVcf.score
        Array[File] prs_raw_scores_log = ScoreVcf.log
        Array[File] prs_sites_scored = ScoreVcf.sites_scored
    }
}