version 1.0

import "../../steps/PRSTasks.wdl"

workflow PRSRawScoreWorkflow {
	input {
		# PRS inputs
		String condition_name
		File var_weights
		File scoring_sites
		File input_vcf
		# Docker images
		String python_docker_image = "python:3.9.10"
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
	}

	call PRSTasks.DetermineChromosomeEncoding as ChrEncoding {
		input:
			weights = var_weights,
			docker_image = python_docker_image
	}

	call PRSTasks.ScoreVcf as ScoreVCF {
		input:
			vcf = input_vcf,
			chromosome_encoding = ChrEncoding.chromosome_encoding,
			sites = scoring_sites,
			weights = var_weights,
			basename = condition_name,
			docker_image = plink_docker_image
	}

	output {
		# Chromosme encoding
		String chromosome_encoding = ChrEncoding.chromosome_encoding
		# PRS outputs
		File prs_raw_scores = ScoreVCF.score
		File prs_raw_scores_log = ScoreVCF.log
    	File prs_sites_scored = ScoreVCF.sites_scored
	}
}