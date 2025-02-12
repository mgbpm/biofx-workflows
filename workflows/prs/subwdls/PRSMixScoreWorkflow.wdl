version 1.0

workflow PRSMixScoreWorkflow {
    input {
        Array[File] raw_scores
        File score_weights
        String output_basename
        # Docker images
        String ubuntu_docker_image = "ubuntu:21.10"
    }

    # Calculate raw PRS mix score
    call CalculateMixScore {
        input:
            raw_scores = raw_scores,
            score_weights = score_weights,
            output_basename = output_basename,
            docker_image = ubuntu_docker_image
    }

    output {
        File prs_mix_raw_score = CalculateMixScore.prs_mix_raw_score
    }
}

task CalculateMixScore {
    input {
        Array[File] raw_scores
        File score_weights
        String output_basename
        String docker_image
        Int disk_size = ceil(size(raw_scores, "GB") + size(score_weights, "GB")) + 10
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        score_file_array=('~{sep="' '" raw_scores}')
        tail -n +2 ${score_file_array[0]} | awk '{ print $1 }' > sample_ids.txt

        printf -- "#IID\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\n" > "~{output_basename}.mix.sscore"

        while read sample_id; do

            score_sum=0
            weight_sum=0

            for score_file in '~{sep="' '" raw_scores}'; do
                pgs_id=$(basename $score_file | grep --ignore-case --only-matching "pgs[0-9]*")
                score_weight=$(grep "${pgs_id}" "~{score_weights}" | cut -f 2)
                raw_score=$(grep "${sample_id}" $score_file | cut -f 4)
                weighted_score=$(awk -v x=${score_weight} -v y=${raw_score} 'BEGIN {print x*y}')
                score_sum=$(awk -v x=${weighted_score} -v y=${score_sum} 'BEGIN {print x+y}')
                weight_sum=$(awk -v x=${score_weight} -v y=${weight_sum} 'BEGIN {print x+y}')
            done

            score_avg=$(awk -v x=${score_sum} -v y=${weight_sum} 'BEGIN {print x/y}')

            printf -- "${sample_id}\t0\t${score_avg}\t${score_sum}\n" >> "~{output_basename}.mix.sscore"

        done < sample_ids.txt
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File prs_mix_raw_score = "~{output_basename}.mix.sscore"
    }
}