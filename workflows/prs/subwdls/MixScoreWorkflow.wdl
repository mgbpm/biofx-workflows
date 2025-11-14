version 1.0

import "../../../steps/Utilities.wdl"

workflow MixScoreWorkflow {
    input {
        Array[File] input_scores
        File        score_weights
        String      output_basename
        String      ubuntu_docker_image = "ubuntu:21.10"
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
        call MixAdjustedScores {
            input:
                input_scores = input_scores,
                score_weights = score_weights,
                output_basename = output_basename,
                docker_image = ubuntu_docker_image
        }
    }
    
    if (is_sscore) {
        call MixRawScores {
            input:
                input_scores = input_scores,
                score_weights = score_weights,
                output_basename = output_basename,
                docker_image = ubuntu_docker_image
        }
    }

    output {
        File mix_score = select_first([MixAdjustedScores.mix_score, MixRawScores.mix_score])
    }
}

task MixRawScores {
    input {
        Array[File] input_scores
        File        score_weights
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
        tail -n+2 ${score_file_array[0]} | awk '{ print $1 }' > sample_ids.txt

        head -n1 "${score_file_array[0]}" > "~{output_basename}.mix.sscore"

        while read sample_id; do
            score_sum=0

            for score_file in '~{sep="' '" input_scores}'; do
                pgs_id=$(basename $score_file | grep --ignore-case --only-matching "pgs[0-9]*")
                score_weight=$(grep "${pgs_id}" "~{score_weights}" | cut -f2)
                raw_score=$(grep "${sample_id}" ${score_file} | cut -f4)
                weighted_score=$(awk -v x=${score_weight} -v y=${raw_score} 'BEGIN {print x*y}')
                score_sum=$(awk -v x=${weighted_score} -v y=${score_sum} 'BEGIN {print x+y}')
            done

            printf -- "${sample_id}\tNA\tNA\t${score_sum}\n" >> "~{output_basename}.mix.sscore"
        done < sample_ids.txt
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{final_disk_size} SSD"
        memory: "~{mem_size} GB"
        preemptible: preemptible
    }

    output {
        File mix_score = "~{output_basename}.mix.sscore"
    }
}

task MixAdjustedScores {
    input {
        Array[File] input_scores
        File        score_weights
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
        tail -n+2 ${score_file_array[0]} | awk '{ print $1 }' > sample_ids.txt

        head -n1 "${score_file_array[0]}" > "~{output_basename}.mix.tsv"

        while read sample_id; do
            score_sum=0

            for score_file in '~{sep="' '" input_scores}'; do
                pgs_id=$(basename ${score_file} | grep --ignore-case --only-matching "pgs[0-9]*")
                score_weight=$(grep "${pgs_id}" "~{score_weights}" | cut -f2)
                adjusted_score=$(grep "${sample_id}" ${score_file} | cut -f26)
                weighted_score=$(awk -v x=${score_weight} -v y=${adjusted_score} 'BEGIN {print x*y}')
                score_sum=$(awk -v x=${weighted_score} -v y=${score_sum} 'BEGIN {print x+y}')
            done

            printf -- "NA\t${sample_id}\t" >> "~{output_basename}.mix.tsv"
            for i in {1..23}; do
                printf -- "NA\t" >> "~{output_basename}.mix.tsv"
            done
            printf -- "${score_sum}\tNA\n" >> "~{output_basename}.mix.tsv"
        done < sample_ids.txt
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{final_disk_size} SSD"
        memory: "~{mem_size} GB"
        preemptible: preemptible
    }

    output {
        File mix_score = "~{output_basename}.mix.tsv"
    }
}