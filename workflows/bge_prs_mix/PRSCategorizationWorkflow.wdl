version 1.0

import "../../steps/PRSTasks.wdl"

workflow PRSCategorizationWorkflow {
	input {
		# PRS inputs
		Array[File] prs_scores
		File sample_ids
		# Ubuntu Docker image
		String ubuntu_docker_image = "ubuntu:21.10"
	}

	call CategorizeScores {
		input:
			prs_scores = prs_scores,
			sample_ids = sample_ids,
			docker_image = ubuntu_docker_image
	}

	output {
		# Individual Outputs
		Array[File] individuals_risk_summaries = CategorizeScores.individual_risk_summaries
	}
}

task CategorizeScores {
	input {
		Array[File] prs_scores
		File sample_ids
		String docker_image
		Int disk_size = ceil(size(prs_scores, "GB")) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		mkdir -p OUTPUT

		while read line; do
			# Write headers for individual's risk summary
			printf "Condition,Percentile,Total_Bins,Individuals_Bin\n">> OUTPUT/"${line}"_risk_summary.csv	

			for c in '~{sep="' '" prs_scores}'; do
				# Find the bin for the percentile
				percentile=$(grep "${line}" $c | cut -f 27)
				if [ $percentile -le 20 ]; then
					bin="Below Average"
				elif [ $percentile -le 75]; then
					bin="Average"
				elif [ $percentile -le 90]; then
					bin="Above Average"
				else
					bin="High"
				fi
				
				# Output individual's percentile and bin to csv
				printf "$(basename ${c} .tsv),${percentile},4,${bin}\n" >> OUTPUT/"${line}"_risk_summary.csv
			done
		done < "~{sample_ids}"
	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		Array[File] individual_risk_summaries = glob("OUTPUT/*.csv")
	}	
}