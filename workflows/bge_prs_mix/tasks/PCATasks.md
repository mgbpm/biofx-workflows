# PCA Tasks for PRS

## MakePCAPlot Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | population_pcs | Yes | Population PCA projections | |
| File | target_pcs | Yes | PCA projections (from flash PCA) | |
| String | output_basename | Yes | Basename for file output | |
| String | docker_image | Yes | Docker image used for running task | "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1" |
| Int | disk_size | No | Disk size to allocate in GB | 100  |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | pca_plot | Always | R-generated plot of input PCs |

## PerformPCA Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | bed | Yes | PLINK's genotype calls | |
| File | bim | Yes | PLINK's extended MAP file with variant information to accompany the input BED file | |
| File | fam | Yes | PLINK's sample information file to accompany the input BED file | |
| String | basename | Yes | Output file base name | |
| String | docker_image | Yes | Docker image used for running task | "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f" |
| Int | nthreads | No | Number of threads to run flash PCA | 16 |
| Int | disk_size | No | Disk size to allocate in GB | 400 |
| Int | mem_size | No | Allocated memory in GB | 8 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | pcs | Always | PCA projections |
| File | pc_variance | Always | PCA variance |
| File | pc_loadings | Always | |
| File | mean_sd | Always | |
| File | eigenvectors | Always | |
| File | eigenvalues | Always | |

## ProjectArray Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | bed | Yes | PLINK's genotype calls | |
| File | bim | Yes | PLINK's extended MAP file with variant information to accompany the input BED file | |
| File | fam | Yes | PLINK's sample information file to accompany the input BED file | |
| File | pc_loadings | Yes | | |
| File | pc_meansd | Yes | | |
| String | basename | Yes | File base name for all output files | |
| Boolean | orderIndependentCheck | No | Post processing step | false |
| String | docker_image | Yes | Docker image used for running task | "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f" |
| Int | nthreads | No | Number of threads to run flash PCA | 16 |
| Int | disk_size | No | Disk size to allocate in GB | 400  |
| Int | mem_size | No | Allocated memory in GB | 8 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | projections | Always | PCA projections output from flash PCA |

## ArrayVcfToPlinkDataset Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | vcf | Yes | VCF to run PLINK2 | |
| File | pruning_sites | Yes | Subset of variants to use for PCA | |
| File | subset_to_sites | No | Subset sites for PCA | |
| String | basename | Yes | File base name for all output files | |
| String | chromosome_encoding | Yes | Chromosome encoding based on if mitochondrial variants are represented by 'chrM' or 'chrMT' | |
| Boolean | use_ref_alt_for_ids | No | | false |
| String | docker_image | Yes | Docker image used for running task | |
| Int | disk_size | No | Disk size to allocate in GB | size of input vcf (x3), plus 20 |
| Int | mem | No | Allocated memory in GB | 8 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | bed | Always | PLINK's genotype calls |
| File | bim | Always | PLINK's extended MAP file with variant information to accompany the output BED file |
| File | fam | Always | PLINK's sample information file to accompany the BED file |
| Array[File] | INPUTS | Always | Summary of input parameters |