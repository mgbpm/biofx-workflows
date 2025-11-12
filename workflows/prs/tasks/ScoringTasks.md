# PRS Scoring Tasks

## ScoreVcf

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

## CheckWeightsCoverSitesUsedInTraining

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | sites_used_in_training | Yes | vars file with sites used in training the population/ancestry adjustment model | |
| WeightSet | weight_set | Yes | Struct of linear weights file, an optional interaction weights file, and SelfExclusiveSites object | |
| String | docker_image | No | Docker image for interaction with python base | "python:3.9.10" |
| Int | disk_size | No | Disk size to allocate in GB | size of sites_used_in_training plus 10 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

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

## AdjustScores

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
