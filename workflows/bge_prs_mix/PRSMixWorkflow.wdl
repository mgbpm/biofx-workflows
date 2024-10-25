version 1.0

import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/ScoringTasks.wdl" as ScoringTasks
import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/PCATasks.wdl" as PCATasks

workflow PRSMixWorkflow {
	input {
		# Scoring inputs
		File input_vcf
		File input_vcf_idx
		Array[File] var_weights
		String population_basename
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
		Array[File] raw_scores
		Int raw_scores_len = length(raw_scores)
		File score_weights
		String output_basename
		String docker_image
		Int disk_size = ceil(size(score_weightss_file, "GB")) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<

		sum_of_score_sum_col=0

		for c in '~{sep="' '" raw_scores}'
		do
			# Get the PGS ID and "IID" column from the var weight file
			pgs_id=$(basename $c .txt | cut -d "_" -f 1)
			iid=$(head -n 2 $c | tail -1 | cut -d "\t" -f 1)
			# Get the weight for the PGS ID in the score weight file
			score_weight=$(grep ${pgs_id} "~{score_weights}" | cut -d "\t" -f 1)
			# Extract the "NAMED_ALLELE_DOSAGE_SUM"/2nd column from the raw score file
			ad_sum=$(head -n 2 $c | tail -1 | cute -d "\t" -f 2)
			# Extract the "SCORE1_AVG"/3rd column from raw score file
			score_avg=$(head -n 2 $c | tail -1 | cute -d "\t" -f 3)
			# Extract the "SCORE1_SUM"/4th column from raw score file
			score_sum=$(head -n 2 $c | tail -1 | cute -d "\t" -f 4)
			# Multiply the PRS raw score by weight and add sum counter
			sum_of_score_sum_col=$((($score_sum * $score_weight) + $sum_of_score_sum_col))
		done

		# Get weighted average for score_sum columns
		weighted_score_sum=$(($sum_of_score_sum_col / "~{raw_scores_len}"))

		# Report the weighted averages in a format similar to raw scores files
		printf "#IID\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\n" > "~{output_basename}_prs_mix_score.sscore"
		printf "${iid}\t${ad_sum}\t${score_avg}\t${weighted_score_sum}\n" >> "~{output_basename}_prs_mix_score.sscore"

	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		File prs_mix_raw_score = "~{output_basename}_prs_mix_score.sscore"
	}
}
