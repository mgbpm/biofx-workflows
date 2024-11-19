version 1.0

import "../../steps/PRSTasks.wdl"

workflow PRSRawScoreWorkflow {
	input {
		# PRS inputs
		Array[File] var_weights
		File scoring_sites
		File input_vcf
		Int? scoring_mem
		# Docker images
		String python_docker_image = "python:3.9.10"
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
	}

	scatter (i in range(length(var_weights))) {
		call PRSTasks.DetermineChromosomeEncoding as ChrEncoding {
			input:
				weights = var_weights[i],
				docker_image = python_docker_image
		}

		call PRSTasks.ScoreVcf as ScoreVCF {
			input:
				vcf = input_vcf,
				chromosome_encoding = ChrEncoding.chromosome_encoding,
				sites = scoring_sites,
				weights = var_weights[i],
				basename = sub(basename(var_weights[i]), ".var_weights.tsv", ""),
				base_mem = select_first([scoring_mem])
				docker_image = plink_docker_image
		}
	}

	output {
		# Chromosme encoding
		Array[String] chromosome_encoding = ChrEncoding.chromosome_encoding
		# PRS outputs
		Array[File] prs_raw_scores = ScoreVCF.score
		Array[File] prs_raw_scores_log = ScoreVCF.log
    	Array[File] prs_sites_scored = ScoreVCF.sites_scored
	}
}