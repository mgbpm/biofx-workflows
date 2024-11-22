version 1.0

workflow PRSMixScoreWorkflow {
	input {
		# PRS Mix inputs
		String condition_name
		Array[File] raw_scores
		File score_weights
		# Docker images
		String ubuntu_docker_image = "ubuntu:21.10"
	}

	# Calculate raw PRS mix score
	call CalculateMixScore {
		input:
			raw_scores = raw_scores,
			score_weights = score_weights,
			output_basename = condition_name,
			docker_image = ubuntu_docker_image
	}

	output {
		File prs_mix_raw_score = CalculateMixScore.prs_mix_raw_score
	}
}

task CalculateMixScore {
	input {
		Array[File] raw_scores
		Int raw_scores_len = length(raw_scores)
		File score_weights
		String output_basename
		String docker_image
		Int disk_size = ceil(size(raw_scores, "GB") + size(score_weights, "GB")) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail


		# Extract all sample IDs from a raw score file
		score_file_array=('~{sep="' '" raw_scores}')
		sed '1d;' ${score_file_array[0]} | awk '{ print $1 }' > sample_ids.txt

		# Set up score file headers
		printf "#IID\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\n" > "~{output_basename}_prs_mix_score.sscore"

		while read line; do
			# Initialize sum of sample's raw scores
			score_sum=0

			# Add the raw score from each file to the sum of raw scores
			for c in '~{sep="' '" raw_scores}'; do
				pgs_id=$(basename $c .var_weights.tsv | cut -d "_" -f 2)
				score_weight=$(grep "${pgs_id}" "~{score_weights}" | cut -f 2)
				raw_score=$(grep "${line}" $c | cut -f 4)
				weighted_score=$(awk -v x=${score_weight} -v y=${raw_score} 'BEGIN {print x*y}')
				score_sum=$(awk -v x=${weighted_score} -v y=${score_sum} 'BEGIN {print x+y}')
			done

			# Get the weighted average for raw scores (the PRS mix score)
			weighted_avg=$(awk -v x=${score_sum} -v y="~{raw_scores_len}" 'BEGIN {print x/y}')

			# Print info for the sample
			printf "${line}\t0\t${weighted_avg}\t${score_sum}\n" >> "~{output_basename}.prs_mix_score.sscore"

		done < sample_ids.txt
	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		File prs_mix_raw_score = "~{output_basename}.prs_mix_score.sscore"
	}
}