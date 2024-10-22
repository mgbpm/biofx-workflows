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
		File score_weight
		# Adjustment inputs
		File population_loadings
		File population_meansd
		File population_pcs
		File pruning_sites_for_pca
		String ubuntu_docker_image = "ubuntu:21.10"
	}

	# Scatter over the files within the zipped file and run scoring on each
	scatter (i in range(length(var_weights))) {
		String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + population_basename
		call ScoringTasks.DetermineChromosomeEncoding as DetermineChrEncoding {
			input:
				weights = var_weights[i]
		}
		call ScoringTasks.ScoreVcf as ScoreImputedVCF {
			input:
				vcf = input_vcf,
				basename = output_basename,
				weights = var_weights[i],
				base_mem = 16,
				chromosome_encoding = DetermineChrEncoding.chromosome_encoding,
		}
	}

	# Calculate raw PRS Mix score
	call CalculateMixScore {
		input:
			raw_scores = ScoreImputedVCF.score,
			score_weight = score_weight
	}

	# Adjust raw PRS Mix score
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
	call ScoringTasks.AdjustScores {
		input:
			fitted_model_params = fitted_params_for_model,
			pcs = ProjectArray.projections,
			scores = CalculateMixScore.prs_mix_raw_score
	}
	call ScoringTasks.MakePCAPlot {
		input:
			population_pcs = population_pcs,
			target_pcs = ProjectArray.projections
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
		File score_weight
		String docker_image
		Int disk_size = ceil(size(score_weights_file, "GB")) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<


	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		File prs_mix_raw_score = "temp"
	}
}
