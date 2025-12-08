version 1.0

import "../../../steps/Utilities.wdl"

workflow MixScoreWorkflow {
    input {
        Array[File] input_scores
        File        score_weights
        String      output_basename
        String      prs_docker_image
    }

    Boolean is_tsv    = basename(input_scores[0]) != basename(input_scores[0], ".tsv")
    Boolean is_sscore = basename(input_scores[0]) != basename(input_scores[0], ".sscore")

    if (!is_tsv && !is_sscore) {
        call Utilities.FailTask {
            input:
                error_message = "Check input score file type. Should be a .tsv or .sscore file."
        }
    }

    if (is_tsv) {
        String tsv_ext = "tsv"
    }
    if (is_sscore) {
        String sscore_ext = "sscore"
    }

    call MixScores {
        input:
            input_scores = input_scores,
            score_weights = score_weights,
            output_ext = select_first([tsv_ext, sscore_ext]),
            output_basename = output_basename,
            docker_image = prs_docker_image
    }

    output {
        File mix_score = MixScores.mix_score
    }
}

task MixScores {
    input {
        Array[File] input_scores
        File        score_weights
        String      output_ext
        String      output_basename
        String      docker_image
        Int         addldisk = 10
        Int         mem_size = 4
        Int         preemptible = 1
    }

    Int file_size       = ceil(size(input_scores, "GB") + size(score_weights, "GB"))
    Int final_disk_size = file_size + addldisk

    command <<<
        set -euxo pipefail

        score_file_array=('~{sep="' '" input_scores}')

        $PACKAGESDIR/biofx-prs/bin/mix_scores.py \
            --score-file-list ${score_file_array} \
            --score-weights "~{score_weights}" \
            --output-basename "~{output_basename}"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{final_disk_size} SSD"
        memory: "~{mem_size} GB"
        preemptible: preemptible
    }

    output {
        File mix_score = "~{output_basename}.mixed.~{output_ext}"
    }
}
