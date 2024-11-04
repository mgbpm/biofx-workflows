# PRS Tasks

## DetermineChromosomeEncoding Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_zip_file | Yes | Tar file of condition-/disease-specific variant weights, score weights, and pca files | |
| String | docker_image | Yes | Docker image used for running task | |
| Int | disk_size | No | Disk size to allocate in GB | size of condition_zip_file plus 10 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| String | chr_encoding | Always | Chromosome encoding based on if mitochondrial variants are represented by 'chrM' or 'chrMT' |

## ScoreVCF Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | VCF to run scoring | |
| String | chr_encoding | Yes | Chromosome encoding based on if mitochondrial variants are represented by 'chrM' or 'chrMT' | |
| File | condition_zip_file | Yes | Tar file of condition-/disease-specific variant weights, score weights, and pca files | |
| String | output_basename | Yes | File base name for all output files | |
| String | docker_image | Yes | Docker image used for running task | |
| Int | disk_size | No | Disk size to allocate in GB | size of input vcf (x3) plus the size of condition_zip_file, plus 20  |
| Int | base_mem | No | Amount of base memory allocation | 16 |
| Int | mem_size | No | Allocated memory in GB | base_mem plus 2 |
| Int | plink_mem | No | Allocated memory in GB for Plink2 | base_mem * 0.75 * 1000 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | score | Always | PRS scores for all samples in the input VCF |
| File | log | Always | Log file for scoring |
| File | sites_scored | Always | Sites from input VCF that were scored |

## ArrayVCFToPlinkDataset Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | VCF to run PLINK2 | |
| File | pruning_sites | Yes | Subset of variants to use for PCA | |
| String | chr_encoding | Yes | Chromosome encoding based on if mitochondrial variants are represented by 'chrM' or 'chrMT' | |
| String | output_basename | Yes | File base name for all output files | |
| String | docker_image | Yes | Docker image used for running task | |
| Int | disk_size | No | Disk size to allocate in GB | size of input vcf (x3), plus 20 |
| Int | mem_size | No | Allocated memory in GB | 8 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | bed | Always | PLINK's genotype calls |
| File | bim | Always | PLINK's extended MAP file with variant information to accompany the output BED file |
| File | fam | Always | PLINK's sample information file to accompany the BED file |

## ProjectArray Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | bed | Yes | PLINK's genotype calls | |
| File | bim | Yes | PLINK's extended MAP file with variant information to accompany the input BED file | |
| File | fam | Yes | PLINK's sample information file to accompany the input BED file | |
| File | condition_zip_file | Yes | Tar file of condition-/disease-specific variant weights, score weights, and pca files | |
| String | output_basename | Yes | File base name for all output files | |
| String | docker_image | Yes | Docker image used for running task | |
| Int | disk_size | No | Disk size to allocate in GB | 400  |
| Int | mem_size | No | Allocated memory in GB | 8 |
| Int | nthreads | No | Number of threads to run flash PCA | 16 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | projections | Always | PCA projections output from flash PCA |

## MakePCAPlot Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_zip_file | Yes | Tar file of condition-/disease-specific variant weights, score weights, and pca files | |
| File | target_pcs | Yes | PCA projections (from flash PCA) | |
| String | output_basename | Yes | File base name for all output files | |
| String | docker_image | Yes | Docker image used for running task | |
| Int | disk_size | No | Disk size to allocate in GB | 100  |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | pca_plot | Always | R-generated plot of input PCs |

## AdjustScores Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_zip_file | Yes | Tar file of condition-/disease-specific variant weights, score weights, and pca files | |
| File | pcs | Yes | PCA projections (from flash PCA) | |
| File | scores | Yes | raw PRS scores | |
| String | output_basename | Yes | File base name for all output files | |
| String | docker_image | Yes | Docker image used for running task | |
| Int | disk_size | No | Disk size to allocate in GB | 100  |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | adjusted_scores | Always | TSV file of PRS scores adjusted using population model |
