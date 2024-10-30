version 1.0

workflow PRSMixWorkflow {
	input {
		# Scoring inputs
		File input_vcf
		File input_vcf_idx
		Array[File] var_weights_files
		# PRS Mix inputs
		File score_weights_file
		# Condition-specific adjustment inputs
		File population_loadings
		File population_meansd
		File population_pcs
		File ancestry_adjustment_model
		# Other adjustment inputs
		File pruning_sites_for_pca
		File scoring_sites
		# Docker images
		String ubuntu_docker_image = "ubuntu:21.10"
	}

	# Calculate raw PRS score per variant weights file
	scatter (i in range(length(var_weights_files))) {
		String var_weight_basename = sub(basename(var_weights_files), ".txt", "")
		call DetermineChromosomeEncoding {
			input:
				weights = var_weights_files[i],
				docker_image = python_docker_image
		}
		call ScoreVcf {
			input:
				input_vcf = input_vcf,
				var_weights = var_weights_files[i],
				chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding,
				sites = scoring_sites,
				output_basename = var_weight_basename + "_" + sample_basename
		}
	}

	# Calculate raw PRS mix score
	call CalculateMixScore {
		input:
			raw_scores = ScoreVcf.score,
			score_weights = score_weights_file,
			output_basename = sample_basename
	}

	# Calculate PCA for individual
	call ExtractIDsPlink {
		input:
			input_vcf = input_vcf,
			chromosome_encoding = DetermineChrEncoding.chromosome_encoding,
			mem = 8
	}
	call ArrayVcfToPlinkDataset {
		input:
		vcf = imputed_array_vcf,
		pruning_sites = pruning_sites_for_pca,
		basename = basename,
		mem = 8,
		use_ref_alt_for_ids = use_ref_alt_for_ids,
		chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding
	}
	call ProjectArray {
		input:
			pc_loadings = population_loadings
			pc_meansd = population_meansd,
			bed = ArrayVcfToPlinkDataset.bed,
			bim = ArrayVcfToPlinkDataset.bim,
			fam = ArrayVcfToPlinkDataset.fam,
			basename = "pca"
	}
	call MakePCAPlot {
		input:
			population_pcs = population_pcs,
			target_pcs = ProjectArray.projections
	}

	# Adjust score with model and PCA
	call AdjustScores {
		input:
			fitted_model_params = ancestry_adjustment_model,
			pcs = ProjectArray.projections,
			scores = CalculateMixScore.prs_mix_raw_score
	}
	
	output {
		# Sample IDs list
		File sample_ids_list = CalculateMixScore.sample_ids_list
		# PRS output
		Array[File] prs_raw_scores = ScoreVcf.score
		Array[File] prs_raw_scores_log = ScoreVcf.log
    	Array[File] prs_raw_sites_scored = ScoreVcf.sites_scored
		File prs_mix_raw_score = CalculateMixScore.prs_mix_raw_score
		File adjusted_scores = AdjustScores.adjusted_scores
		# PCA output
		File pc_projection = ProjectArray.projections
		File pc_plot = MakePCAPlot.pca_plot
	}
}

task CalculateMixScore {
	input {
		Array[File] raw_scores_files
		Int raw_scores_len = length(raw_scores_files)
		File score_weights
		String output_basename
		String docker_image
		Int disk_size = ceil(size(score_weights, "GB") * 2) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		mkdir -p OUTPUT
		mkdir -p OUTPUT

		# Extract all sample IDs from a raw score file
		score_file_array=('~{sep="' '" raw_scores_files}')
		sed '1d;' ${score_file_array[0]} | awk '{ print $1 }' > OUTPUT/sample_ids.txt

		# Set up score file headers
		printf "#IID\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\n" > "OUTPUT/~{output_basename}_prs_mix_score.sscore"

		while read line; do
			# Initialize sum of sample's raw scores
			score_sum=0

			# Add the raw score from each file to the sum of raw scores
			for c in '~{sep="' '" raw_scores_files}'; do
				pgs_id=$(basename $c .txt | cut -d "_" -f 1)
				score_weight=$(grep "${pgs_id}" "~{score_weights}" | cut -f 2)
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

task ScoreVcf {
	input {
    	File vcf
    	String basename
    	File weights
    	Int base_mem = 8
    	String? extra_args
    	File? sites
    	File? exclude_sites
    	String? chromosome_encoding
    	Boolean use_ref_alt_for_ids = false
    	Boolean use_dosage_annotation = true
  	}

  	Int runtime_mem = base_mem + 2
  	Int plink_mem = ceil(base_mem * 0.75 * 1000)
  	Int disk_space =  3*ceil(size(vcf, "GB")) + 20
  	String var_ids_string = "@:#:" + if use_ref_alt_for_ids then "\\$r:\\$a" else "\\$1:\\$2"

  	command <<<
  		/plink2 --score ~{weights} header ignore-dup-ids list-variants no-mean-imputation \
  		cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums --set-all-var-ids ~{var_ids_string} --allow-extra-chr ~{extra_args} -vcf ~{vcf} ~{if use_dosage_annotation then "dosage=DS" else ""} \
  		--new-id-max-allele-len 1000 missing ~{"--extract " + sites} ~{"--exclude " + exclude_sites} --out ~{basename} --memory ~{plink_mem} ~{"--output-chr " + chromosome_encoding}
  	>>>

  	output {
  	  	File score = "~{basename}.sscore"
  	  	File log = "~{basename}.log"
  	  	File sites_scored = "~{basename}.sscore.vars"
  	}

  	runtime {
  		docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
  		disks: "local-disk " + disk_space + " HDD"
  		memory: runtime_mem + " GB"
  	}
}

task ExtractIDsPlink {
	input {
		File vcf
		Int disk_size = 2 * ceil(size(vcf, "GB")) + 100
		Boolean use_ref_alt_for_ids = false
		Int mem = 8
		String? chromosome_encoding
	}

  	Int plink_mem = ceil(mem * 0.75 * 1000)
  	String var_ids_string = "@:#:" + if use_ref_alt_for_ids then "\\$r:\\$a" else "\\$1:\\$2"

  	command <<<
  		/plink2 \
  			--vcf ~{vcf} \
  			--set-all-var-ids ~{var_ids_string} \
  			--new-id-max-allele-len 1000 missing \
  			--rm-dup exclude-all \
  			--allow-extra-chr \
  			--write-snplist allow-dups \
  			--memory ~{plink_mem} ~{"--output-chr " + chromosome_encoding}
  	>>>

	output {
	  	File ids = "plink2.snplist"
	}

  	runtime {
  	  	docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
  	  	disks: "local-disk " + disk_size + " HDD"
  	  	memory: mem + " GB"
  	}
}

#plink chromosome encoding rules: https://www.cog-genomics.org/plink/2.0/data#irreg_output
task DetermineChromosomeEncoding {
  	input {
  	  	File weights
  	}

  	command <<<
  	  	python3 << "EOF"
  	  	code = 'MT'
  	  	with open("~{weights}") as weights_file:
  	  	  	chroms = {s.split("\t")[0].split(":")[0] for i, s in enumerate(weights_file) if i > 0}
  	  	  	if any('chr' in c for c in chroms):
  	  	  	    if 'chrM' in chroms:
  	  	  	        code = 'chrM'
  	  	  	    else:
  	  	  	        code = 'chrMT'
  	  	  	elif 'M' in chroms:
  	  	      	code = 'M'

  	  	with open('chr_encode_out.txt', 'w') as write_code_file:
  	  		write_code_file.write(f'{code}\n')
  	  	EOF
  	>>>

  	runtime {
  	  docker : "python:3.9.10"
  	}

  	output {
  	  String chromosome_encoding = read_string("chr_encode_out.txt")
  }
}

task ArrayVcfToPlinkDataset {
	input {
	  File vcf
	  File pruning_sites
	  File? subset_to_sites
	  String basename
	  Int mem = 8
	  Boolean use_ref_alt_for_ids = false
	  String? chromosome_encoding
	}

	Int disk_space =  3 * ceil(size(vcf, "GB")) + 20
	String var_ids_string = "@:#:" + if use_ref_alt_for_ids then "\\$r:\\$a" else "\\$1:\\$2"	
	
	command <<<
	  	/plink2 --vcf ~{vcf} --extract-intersect ~{pruning_sites} ~{subset_to_sites} --allow-extra-chr --set-all-var-ids ~{var_ids_string} \
	  	--new-id-max-allele-len 1000 missing --out ~{basename} --make-bed --rm-dup force-first ~{"--output-chr " + chromosome_encoding}
	>>>

	output {
	  	File bed = "~{basename}.bed"
	  	File bim = "~{basename}.bim"
	  	File fam = "~{basename}.fam"
	}	
	runtime {
	  	docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
	  	disks: "local-disk " + disk_space + " HDD"
		memory: mem + " GB"
	}
}

task ProjectArray {
  	input {
  	  	File bim
  	  	File bed
  	  	File fam
  	  	File pc_loadings
  	  	File pc_meansd
  	  	String basename
  	  	Int mem = 8
  	  	Int nthreads = 16
  	  	String? divisor
  	}

  	command <<<
  	  	cp ~{bim} ~{basename}.bim
  	  	cp ~{bed} ~{basename}.bed
  	  	cp ~{fam} ~{basename}.fam

  	  	cp ~{pc_loadings} loadings.txt
  	  	cp ~{pc_meansd} meansd.txt

  	  	# Check if .bim file, pc loadings, and pc meansd files have the same IDs
  	  	# 1. extract IDs, removing first column of .bim file and first rows of the pc files
  	  	awk '{print $2}' ~{basename}.bim > bim_ids.txt
  	  	awk '{print $1}' loadings.txt | tail -n +2 > pcloadings_ids.txt
  	  	awk '{print $1}' meansd.txt | tail -n +2 > meansd_ids.txt

  	  	diff bim_ids.txt pcloadings_ids.txt > diff1.txt
  	  	diff bim_ids.txt meansd_ids.txt > diff2.txt
  	  	diff pcloadings_ids.txt meansd_ids.txt > diff3.txt

  	  	if [[ -s diff3.txt ]]
  	  	then
  	  	echo "PC loadings file and PC means file do not contain the same IDs; check your input files and run again."
  	  	exit 1
  	  	fi

  	  	# check if diff files are not empty
  	  	if [[ -s diff1.txt || -s diff2.txt ]]
  	  	then
  	  	echo "IDs in .bim file are not the same as the IDs in the PCA files; check that you have the right files and run again."
  	  	exit 1
  	  	fi

  	  	~/flashpca/flashpca \
  	  	  --bfile ~{basename} \
  	  	  --numthreads ~{nthreads} \
  	  	  --project \
  	  	  --inmeansd meansd.txt \
  	  	  --outproj projections.txt \
  	  	  --inload loadings.txt \
  	  	  -v \
  	  	  ~{"--div " + divisor}
  	>>>

  	output {
  	  	File projections = "projections.txt"
  	}

  	runtime {
  	  	docker: "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
  	  	disks: "local-disk 400 HDD"
  	  	memory: mem + " GB"
  	}
}

task MakePCAPlot {
  	input {
  	  	File population_pcs
  	  	File target_pcs
  	}

  	command <<<
  	  	Rscript -<< "EOF"
  	  	  	library(dplyr)
  	  	  	library(readr)
  	  	  	library(ggplot2)

  	  	  	population_pcs <- read_tsv("~{population_pcs}")
  	  	  	target_pcs <- read_tsv("~{target_pcs}")

  	  	  	ggplot(population_pcs, aes(x=PC1, y=PC2, color="Population")) +
  	  	  		geom_point(size=0.1, alpha=0.1) +
  	  	  	  	geom_point(data = target_pcs, aes(x=PC1, y=PC2, color="Target")) +
  	  	  	  	labs(x="PC1", y="PC2") +
  	  	  	  	theme_bw()

  	  	  	ggsave(filename = "PCA_plot.png", dpi=300, width = 6, height = 6)

  	  	EOF
  	>>>

  	output {
  		File pca_plot = "PCA_plot.png"
  	}

  	runtime {
  	  	docker: "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
  	  	disks: "local-disk 100 HDD"
  	}
}

task AdjustScores {
  	input {
		File fitted_model_params
		File pcs
		File scores
		Int mem = 2
  	}

  	command <<<
  	  Rscript -<< "EOF"
  	    library(dplyr)
  	    library(readr)

  	    # read in model params
  	    params_tibble <- read_tsv("~{fitted_model_params}")
  	    params <- params_tibble %>% pull(value)

  	    # linear transformation to predict variance
  	    f_sigma2 <- function(t, theta) {
  	      PC1 = t %>% pull(PC1)
  	      PC2 = t %>% pull(PC2)
  	      PC3 = t %>% pull(PC3)
  	      PC4 = t %>% pull(PC4)
  	      PC5 = t %>% pull(PC5)
  	      sigma2 <- exp(theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4)
  	    }


  	    # linear transformation to predict mean
  	    f_mu <- function(t, theta) {
  	      PC1 = t %>% pull(PC1)
  	      PC2 = t %>% pull(PC2)
  	      PC3 = t %>% pull(PC3)
  	      PC4 = t %>% pull(PC4)
  	      PC5 = t %>% pull(PC5)
  	      mu <- theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4
  	    }

  	    scores = inner_join(read_tsv("~{pcs}"),
  	                              read_tsv("~{scores}"), by=c("IID" = "#IID"))

  	    adjusted_scores <- scores %>% mutate(adjusted_score =
  	                                                      (SCORE1_SUM - f_mu(scores, params[1:5]))/
  	                                                      sqrt(f_sigma2(scores, params[6:10]))
  	                                                    )
  	    adjusted_scores <- adjusted_scores %>% mutate(percentile=pnorm(adjusted_score,0))

  	    # return array scores
  	    write_tsv(adjusted_scores, "adjusted_scores.tsv")
  	  EOF
  	>>>

  	output {
		File adjusted_scores = "adjusted_scores.tsv"
  	}

  	runtime {
		docker: "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
		disks: "local-disk 100 HDD"
		memory: mem + " GB"
  	}
}