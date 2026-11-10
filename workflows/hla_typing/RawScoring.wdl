version 1.0

import "./ScoringTasks.wdl"

workflow RawScoreWorkflow {
    input {
        File?  input_vcf
        File?  input_pgen
        File?  input_pvar
        File?  input_psam
        File   weights
        String plink2_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    }

    call ScoringTasks.RawPrsTask {
        input:
            vcf     = input_vcf,
            pgen    = input_pgen,
            pvar    = input_pvar,
            psam    = input_psam,
            weights = weights,
            docker_image = plink2_docker_image
    }

    output {
        File prs_score        = RawPrsTask.score
        File sites_scored     = RawPrsTask.sites_scored
        File sites_not_scored = RawPrsTask.sites_not_scored
    }
}
