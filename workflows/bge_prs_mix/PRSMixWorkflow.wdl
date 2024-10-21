version 1.0

import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/ScoringTasks.wdl" as ScoringTasks
import "https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/PCATasks.wdl" as PCATasks

workflow PRSMixWorkflow {
	input {
		# Scoring inputs
		File imputed_vcf
		File imputed_vcf_idx
		Array[File] var_weights
		# PRS Mix inputs
		File score_weight
		# Adjustment inputs
		File population_vcf
		File pruning_sites_for_pca
		String ubuntu_docker_image = "ubuntu:21.10"
	}

	# Scatter over the files within the zipped file and run scoring on each
	scatter (i in range(length(var_weights))) {
		String output_basename = sub(basename(var_weights[i]), "txt", "")
		call ScoringTasks.DetermineChromosomeEncoding as DetermineChrEncoding {
			input:
				weights = var_weights[i]
		}
		call ScoringTasks.ScoreVcf as ScoreImputedVCF {
			input:
				vcf = imputed_vcf,
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
			vcf = imputed_vcf,
			chromosome_encoding = DetermineChrEncoding.chromosome_encoding,
			mem = 8
	}
	call PCATasks.ArrayVcfToPlinkDataset as PopulationArrayVcfToPlinkDataset {
		input:
			vcf = population_vcf,
			pruning_sites = pruning_sites_for_pca,
			subset_to_sites = ExtractIDsPlink.ids,
			basename = "population",
			use_ref_alt_for_ids = use_ref_alt_for_ids,
			chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding
	}
	call PCATasks.PerformPCA {
	  input:
		bim = PopulationArrayVcfToPlinkDataset.bim,
		bed = PopulationArrayVcfToPlinkDataset.bed,
		fam = PopulationArrayVcfToPlinkDataset.fam,
		basename = basename
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
	call CheckPopulationIds {
		input:
			pop_vcf_ids = ExtractIDsPopulation.ids,
			pop_pc_loadings = PerformPCA.pc_loadings,
			docker_image = ubuntu_docker_image
	}
	call PCATasks.ProjectArray {
		input:
			pc_loadings = PerformPCA.pc_loadings
			pc_meansd = PerformPCA.mean_sd,
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
			population_pcs = PerformPCA.pcs,
			target_pcs = ProjectArray.projections
	}
	if (!select_first([CheckPopulationIdsValid.files_are_valid, true])) {
		call ErrorWithMessage {
			input:
			message = "Population VCF IDs are not a subset of the population PCA IDs; running with these inputs would give an incorrect result."
		}
	}

	output {
		# Scores
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

task CheckPopulationIds {
	input {
		File pop_vcf_ids
		File pop_pc_loadings
		String docker_image
	}
	command <<<
		# check if population VCF file contains a subset of population PC loading ids

		# 1. extract IDs, removing first rows of the pc files
		awk '{print $1}' ~{pop_pc_loadings} | tail -n +2 > pop__pc_ids.txt

		comm -23 <(sort pop_pc_ids.txt | uniq) <(sort ~{pop_vcf_ids} | uniq) > array_specific_ids.txt
		if [[ -s array_specific_ids.txt ]]
		then
		echo false
		else
		echo true
		fi

	>>>
	output {
		Boolean files_are_valid = read_boolean(stdout())
	}
	runtime {
		docker: "~{docker_image}"
	}
}