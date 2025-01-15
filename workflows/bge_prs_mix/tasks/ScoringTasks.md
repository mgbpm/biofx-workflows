# PRS Tasks

## ScoreVcf Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | vcf | Yes | VCF to run scoring | |
| String | basename | Yes | Output files base name | |
| File | weights | Yes | Variant weights file | |
| String | extra_args | No | | |
| File | sites | No | | |
| File | exclude_sites | No | | |
| String | chromosome_encoding | No | Chromosome encoding based on if mitochondrial variants are represented by 'chrM' or 'chrMT' | |
| Boolean | use_ref_alt_for_ids | No | | false |
| Boolean | use_dosage_annotation | No | | false |
| String | docker_image | Yes | Docker image used for running task | "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124" |
| Int | base_mem | No | Amount of base memory allocation | 8 |
| Int | disk_size | No | Disk size to allocate in GB | size of input vcf x3, plus 20  |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | score | Always | PRS scores for all samples in the input VCF |
| File | log | Always | Log file for scoring |
| File | sites_scored | Always | Sites from input VCF that were scored |
| Array[File] | INPUTS | Always | Summary of inputs |

## AddInteractionTermsToScore Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | vcf | Yes | VCF with all samples | |
| File | interaction_weights | Yes | | |
| File | scores | Yes | | |
| File | sites | No | File with columns ID, chrom, and pos | |
| String | basename | Yes | Output file base name | |
| SelfExclusiveSites | self_exclusive_sites | No | Struct of sites file (with ID, chrom, and pos columns) and the number of maximum allowed sites | |
| String | docker_image | No | Docker image for interaction with python base | "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e" |
| Int | disk_size | No | Disk size to allocate in GB | 3x the input VCF, plus 20 |
| Int | mem_size | No | Allocated memory in GB | 8 |
| Int | block_buffer | No | | 10000000 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | scores_with_interactions | Always | |
| File | sites_used_in_interaction_score | Always | |

## CheckWeightsCoverSitesUsedInTraining Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | sites_used_in_training | Yes | vars file with sites used in training the population/ancestry adjustment model | |
| WeightSet | weight_set | Yes | Struct of linear weights file, an optional interaction weights file, and SelfExclusiveSites object | |
| String | docker_image | No | Docker image for interaction with python base | "python:3.9.10" |
| Int | disk_size | No | Disk size to allocate in GB | size of sites_used_in_training plus 10 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

## CompareScoredSitesToSitesUsedInTraining Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | sites_used_in_training | Yes | vars file with sites used in training the population/ancestry adjustment model | |
| File | sites_used_in_scoring | Yes | vars file with sites used in finding the raw PRS scores for samples | |
| WeightSet | weight_set | Yes | Struct of linear weights file, an optional interaction weights file, and SelfExclusiveSites object | |
| String | docker_image | No | Docker image for interaction with python base | "python:3.9.10" |
| Int | disk_size | No | Disk size to allocate in GB | size of sites_used_in_training plus size of sites_used_in_scoring, plus 10 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | missing_sites | Always | Sites used in training the model, but not in scoring |
| File | n_missing_sites | Always | Number of missing sites |
| Float | max_error_up | Always |  |
| Float | max_error_down | Always |  |

## CombineScoringSites Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | sites_used_linear_score | Yes | | |
| File | sites_used_interaction_score | Yes | | |
| String | basename | Yes | Output file base name | |
| String | docker_image | No | Ubuntu Docker image | "ubuntu:20.04" |
| Int | disk_size | No | Disk size to allocate in GB | size of sites_used_linear_score plus 2x sites_used_interaction_score, plus 50 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | combined_scoring_sites | Always | Sites from both linear scores and interaction scores files |

## AddShiftToRawScores Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | raw_scores | Yes | Raw PRS scores | |
| Float | shift | Yes | Amount of shift to add to raw PRS scores | |
| String | basename | Yes | Output file base name | |
| String | docker_image | No | tidyverse Docker image | "rocker/tidyverse:4.1.0" |
| Int | disk_size | No | Disk size to allocate in GB | 100 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | shifted_scores | Always | TSV of shifted raw PRS scores |

## CombineMissingSitesAdjustedScores Tasl

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | adjusted_scores_shifted_up | Yes | Shifted up adjusted scores | |
| File | adjusted_scores_shifted_down | Yes | Shifted down adjusted scores | |
| File | adjusted_scores | Yes | PRS raw scores adjusted using population model | |
| Int | n_missing_sites | Yes | Number of sites missing | |
| String | condition_name | Yes | Name of disease/condition | |
| String | docker_image | No | tidyverse Docker image | "rocker/tidyverse:4.1.0" |
| Int | disk_size | No | Disk size to allocate in GB | 100 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | missing_sites_shifted_scores | Always | Shifted scores for the missing sites |

## TrainAncestryModel

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | population_pcs | Yes | Population PCA projections | |
| File | population_scores | Yes | Score for population VCF | |
| String | output_basename | Yes | Output file base name | |
| String | docker_image | No | tidyverse Docker image | "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1" |
| Int | disk_size | No | Disk size to allocate in GB | 100 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | fitted_params | Always | Adjustment model parameters fitted to input population |
| File | adjusted_population_scores | Always | Population scores adjusted with model parameters |
| Boolean | fit_converged | Always |

## AdjustScores Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | fitted_model_params | Yes | Trained adjustment model parameters | |
| File | pcs | Yes | PCA projections (from flash PCA) | |
| File | scores | Yes | Raw PRS scores | |
| String | output_basename | Yes | Basename for file output | |
| String | docker_image | Yes | Docker image used for running task | "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1" |
| Int | disk_size | No | Disk size to allocate in GB | 100  |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | adjusted_scores | Always | TSV file of PRS scores adjusted using population model |

## ExtractIDsPlink

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | vcf | Yes | VCF with all samples | |
| Boolean | use_ref_alt_for_ids | No | | false |
| String | chromosome_encoding | No | Chromosome encoding based on if mitochondrial variants are represented by 'chrM' or 'chrMT' | |
| String | docker_image | Yes | Docker image used for running task | "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124" |
| Int | disk_size | No | Disk size to allocate in GB | 2x the input VCF, plus 100 |
| Int | mem_ | No | Allocated memory in GB | 8 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | ids | Always | Plink2 SNP list |

## DetermineChromosomeEncoding Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | weights | Yes | Variant weights file | |
| String | docker_image | Yes | Docker image used for running task | "python:3.9.10" |
| Int | disk_size | No | Disk size to allocate in GB | size of weights file plus 10 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| String | chromosome_encoding | Always | Chromosome encoding based on if mitochondrial variants are represented by 'chrM' or 'chrMT' |
