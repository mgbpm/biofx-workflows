version 1.0

import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/ScoringTasks.wdl" as ScoringTasks
import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/PCATasks.wdl" as PCATasks

workflow PRSMixWorkflow {
	input {
		# Scoring inputs
		File input_vcf
		File input_vcf_idx
		Array[File] var_weights
		# PRS Mix inputs
		File score_weights
		# Adjustment inputs
		File population_loadings
		File population_meansd
		File population_pcs
		File pruning_sites_for_pca
		String ubuntu_docker_image = "ubuntu:21.10"
	}

	# Scatter over the files within the zipped file and run scoring on each
	String sample_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
	scatter (i in range(length(var_weights))) {
		String var_weight_basename = sub(basename(var_weights), ".txt", "")
		call ScoringTasks.DetermineChromosomeEncoding as DetermineChrEncoding {
			input:
				weights = var_weights[i]
		}
		call ScoringTasks.ScoreVcf as ScoreImputedVCF {
			input:
				vcf = input_vcf,
				basename = var_weight_basename + "_" + sample_basename,
				weights = var_weights[i],
				base_mem = 16,
				chromosome_encoding = DetermineChrEncoding.chromosome_encoding,
		}
	}

	# Calculate raw PRS Mix score
	call CalculateMixScore {
		input:
			raw_scores = ScoreImputedVCF.score,
			score_weights = score_weights,
			output_basename = sample_basename
	}

	# Calculate PCA for individual
	call ScoringTasks.ExtractIDsPlink {
		input:
			vcf = input_vcf,
			chromosome_encoding = DetermineChrEncoding.chromosome_encoding,
			mem = 8
	}
	call PCATasks.ArrayVcfToPlinkDataset {
		input:
		vcf = imputed_array_vcf,
		pruning_sites = pruning_sites_for_pca,
		basename = basename,
		mem = 8,
		use_ref_alt_for_ids = use_ref_alt_for_ids,
		chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding
	}
	call PCATasks.ProjectArray {
		input:
			pc_loadings = population_loadings
			pc_meansd = population_meansd,
			bed = ArrayVcfToPlinkDataset.bed,
			bim = ArrayVcfToPlinkDataset.bim,
			fam = ArrayVcfToPlinkDataset.fam,
			basename = "pca"
	}
	call ScoringTasks.MakePCAPlot {
		input:
			population_pcs = population_pcs,
			target_pcs = ProjectArray.projections
	}

	# Adjust score with model and PCA
	call ScoringTasks.AdjustScores {
		input:
			fitted_model_params = fitted_params_for_model,
			pcs = ProjectArray.projections,
			scores = CalculateMixScore.prs_mix_raw_score
	}
	
	output {
		# PRS Scores
		Array[File] prs_raw_scores = ScoreImputedVCF.score
		Array[File] prs_raw_scores_log = ScoreImputedVCF.log
    	Array[File] prs_raw_sites_scored = ScoreImputedVCF.sites_scored
		File prs_mix_raw_score = CalculateMixScore.prs_mix_raw_score
		File adjusted_scores = AdjustScores.adjusted_scores
		# PCA
		File pc_projection = ProjectArray.projections
		File pc_plot = MakePCAPlot.pca_plot
	}
}

task CalculateMixScore {
	input {
		Array[File] raw_scores_files
		Int raw_scores_len = length(raw_scores_files)
		File score_weights_file
		String output_basename
		String docker_image
		Int disk_size = ceil(size(score_weights_file, "GB") * 2) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		mkdir -p OUTPUT
		mkdir -p WORK

		# Extract all sample IDs from a raw score file
		score_file_array=('~{sep="' '" raw_scores_files}')
		sed '1d;' ${score_file_array[0]} | awk '{ print $1 }' > WORK/sample_ids.txt

		# Set up score file headers
		printf "#IID\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\n" > "OUTPUT/~{output_basename}_prs_mix_score.sscore"

		while read line; do
			# Initialize sum of sample's raw scores
			score_sum=0

			# Add the raw score from each file to the sum of raw scores
			for c in '~{sep="' '" raw_scores_files}'; do
				pgs_id=$(basename $c .txt | cut -d "_" -f 1)
				score_weight=$(grep "${pgs_id}" "~{score_weights_file}" | cut -f 2)
				raw_score=$(grep "${line}" $c | cut -f 4)
				weighted_score=$(awk -v x=${score_weight} -v y=${raw_score} 'BEGIN {print x*y}')
				score_sum=$(awk -v x=${weighted_score} -v y=${score_sum} 'BEGIN {print x+y}')
			done

			# Get the weighted average for raw scores (the PRS mix score)
			weighted_avg=$(awk -v x=${score_sum} -v y="~{raw_scores_len}" 'BEGIN {print x/y}')

			# Print info for the sample
			printf "${line}\t0\t${weighted_avg}\t${score_sum}\n" >> "OUTPUT/~{output_basename}_prs_mix_score.sscore"

		done < WORK/sample_ids.txt
	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		File prs_mix_raw_score = "OUTPUT/~{output_basename}_prs_mix_score.sscore"
	}
}