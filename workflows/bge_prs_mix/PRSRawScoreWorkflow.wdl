version 1.0

import "PRSTasks.wdl"

workflow PRSRawScoreWorkflow {
	input {
		# Input VCF
		File input_vcf
		# Zip with condition-specific files
		File condition_file
		# Docker images
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		String interaction_docker_image = "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e"
	}

	String condition_name = sub(basename(condition_file), "\\.(tar|TAR|tar.gz|TAR.GZ)$", "")

	call PRSTasks.DetermineChromosomeEncoding as ChrEncoding {
		input:
			condition_zip_file = condition_file,
			docker_image = interaction_docker_image
	}
	call PRSTasks.ScoreVCF as GetRawScores {
		input:
			input_vcf = input_vcf,
			chromosome_encoding = ChrEncoding.chr_encoding,
			condition_zip_file = condition_file,
			output_basename = condition_name,
			docker_image = plink_docker_image
	}

	output {
		File chromosome_encoding = ChrEncoding.chr_encoding
		File prs_raw_scores = GetRawScores.score
		File prs_raw_scores_log = GetRawScores.log
    	File prs_raw_sites_scored = GetRawScores.sites_scored
	}
}