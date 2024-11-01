version 1.0

workflow PRSMixScoreWorkflow {
	input {
		# Files with raw scores
		Array[File] raw_scores
		# Zip with condition-specific files
		File condition_file
		# Docker images
		String ubuntu_docker_image = "ubuntu:21.10"
	}

	String condition_name = sub(basename(condition_file), "\\.(tar|TAR|tar.gz|TAR.GZ)$", "")

	# Calculate raw PRS mix score
	call CalculateMixScore {
		input:
			raw_scores = raw_scores,
			condition_zip_file = condition_file,
			output_basename = condition_name,
			docker_image = ubuntu_docker_image
	}

	output {
		File prs_mix_raw_score = CalculateMixScore.prs_mix_raw_score
		File sample_ids_list = CalculateMixScore.sample_ids_list
	}
}

task CalculateMixScore {
	input {
		Array[File] raw_scores
		Int raw_scores_len = length(raw_scores)
		File condition_zip_file
		String output_basename
		String docker_image
		Int disk_size = ceil(size(raw_scores, "GB") + size(condition_zip_file, "GB")) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		mkdir -p OUTPUT
		mkdir -p WORK

		score_weights_file=$(tar -tf "~{condition_zip_file}" | grep score_weights | cut -d "/" -f 2)
		tar -xf "~{condition_zip_file}" ${score_weights_file}
		mv ${score_weights_file} WORK/score_weights.txt

		# Extract all sample IDs from a raw score file
		score_file_array=('~{sep="' '" raw_scores}')
		sed '1d;' ${score_file_array[0]} | awk '{ print $1 }' > OUTPUT/sample_ids.txt

		# Set up score file headers
		printf "#IID\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\n" > "OUTPUT/~{output_basename}_prs_mix_score.sscore"

		while read line; do
			# Initialize sum of sample's raw scores
			score_sum=0

			# Add the raw score from each file to the sum of raw scores
			for c in '~{sep="' '" raw_scores}'; do
				pgs_id=$(basename $c .txt | cut -d "_" -f 1)
				score_weight=$(grep "${pgs_id}" WORK/score_weights.txt | cut -f 2)
				raw_score=$(grep "${line}" $c | cut -f 4)
				weighted_score=$(awk -v x=${score_weight} -v y=${raw_score} 'BEGIN {print x*y}')
				score_sum=$(awk -v x=${weighted_score} -v y=${score_sum} 'BEGIN {print x+y}')
			done

			# Get the weighted average for raw scores (the PRS mix score)
			weighted_avg=$(awk -v x=${score_sum} -v y="~{raw_scores_len}" 'BEGIN {print x/y}')

			# Print info for the sample
			printf "${line}\t0\t${weighted_avg}\t${score_sum}\n" >> "OUTPUT/~{output_basename}_prs_mix_score.sscore"

		done < OUTPUT/sample_ids.txt
	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		File sample_ids_list = "OUTPUT/sample_ids.txt"
		File prs_mix_raw_score = "OUTPUT/~{output_basename}_prs_mix_score.sscore"
	}
}