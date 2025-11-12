# Helper Tasks for PRS

## GetBaseMemory

Computes the memory (in gigabytes) required by plink to process vcf, according to the recommendations given in <https://www.cog-genomics.org/plink/2.0/other#memory>

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---  | :--- | :--- |
| File | vcf | | VCF that will be used in PLINK | |
| Int | nvariants | | Number of variants in VCF | |
| String | docker_image | No | Docker image for VM | "python:3.11" |
| Int | addldisk | No | Extra disk size to allocate in GB | 20 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

**NB** Exactly one of  `vcf` and  `nvariants` must be specified.

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Int | gigabytes | Always | Memory to allocate to PLINK |
| Int | nvariants_ | Always | Number of variants in input VCF |

## RenameChromosomesInTsv

Renames the human chromosomes mentioned in the first column of  `tsv` to 1, 2, ..., 22, X, Y, MT.

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| File | tsv | Yes | TSV with chromosomes to rename | |
| Boolean | skipheader | Yes | Whether or not to skip the header of the TSV | |
| File | lookup | No | Chromosome name mapping | "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv" |
| String | docker_image | No | Docker image for VM | "python:3.11" |
| Int | addldisk | No | Extra disk size to allocate in GB | 20 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | renamed | Always | TSV with renamed chromosomes |

## RenameChromosomesInVcf

Renames the human chromosomes mentioned in the `CHROM` column of `vcf` to 1, 2, ..., 22, X, Y, MT.

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| File | vcf | Yes | VCF with chromosomes to rename | |
| File | lookup | No | Chromosome name mapping | "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv" |
| String | docker_image | No | Docker image for VM | "biocontainers/bcftools:v1.9-1-deb_cv1" |
| Int | addldisk | No | Extra disk size to allocate in GB | 20 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name    | When   | Description    |
| :--- | :---    | :---   | :---    |
| File | renamed | Always | VCF with renamed chromosomes |

## SubsetVcf

Subset  `inputvcf ` to the regions specified in the  `regions` file.

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| File | inputvcf| Yes | VCF for subsetting | |
| File | regions | Yes | Regions for subsetting the input VCF | |
| String | label | No | Basename for the subset VCF | "data" |
| Boolean | nocleanup | No | Include clean-up steps in the process | false |
| String | docker_image | No | Docker image for VM | "biocontainers/bcftools:v1.9-1-deb_cv1" |
| Int | addldisk | No | Extra disk size to allocate in GB | 20 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | result | Always | Subsetted VCF |
| Int | nregions | Always | Number of regions in VCF |

## Union

Generate the set union of the lines in all the files in  `lists`.  The task simply returns the output of

    sort --unique FILE_0 FILE_1 ... FILE_n

...where  `FILE_0`,  `FILE_1`, ...,  `FILE_n` are the files in  `lists`

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---  | :--- | :--- |
| Array[File]+ | lists | Yes | (Nonempty) array of text files | |
| String | docker_image | No | Docker image for VM | "ubuntu:21.10" |
| Int | addldisk | No | Extra disk size to allocate in GB | 20 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | result | Always | File holding the union of lines in  `lists` |

## Intersection

Generate the set intersection of the lines in all the files in
 `lists `.

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| Array[File]+ | lists | Yes | (Nonempty) array of text files | |
| String | docker_image | No | Docker image for VM | "ubuntu:21.10" |
| Int | addldisk | No | Extra disk size to allocate in GB | 20 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | result | Always | File holding the intersection of lines  `lists` |

## ListShards

List all the VCF shard files in the directory tree rooted at
 `source `.  (Here, "VCF shard files" means files with
extension ".vcf.gz" whose parent directory's basename is "shards".)

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---  | :--- | :--- |
| String | source | Yes | URL to a GCP bucket | |
| String | workspace | Yes | Name of Terra workspace | |
| String | docker_image | No | Docker image for VM | |
| Int | disk_size | No | Disk size to allocate in GB | 10 |
| Int | mem_size | No | Allocated memory in GB | 2 |
| Int | preemptible | No | Preemptible runtime setting | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :---   | :--- |
| File | relpaths | Always | File the paths of all VCF shard paths relative to  `source` |

## MakeBatches

Distribute the items in  `cases ` across
 `nbatches ` batches of roughly equal size.

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---  | :--- | :--- |
| File | cases | Yes | File containing cases to batch, one per line | |
| Int | nbatches | Yes | Desired number of batches | |
| Array[String] | exclude | No | (Possibly empty) array of exclusions | [] |
| Boolean| noshuffle | No | If true, do not shuffle the cases before batching them | false |
| Int? | seed | No | Seed to use for the PRNG before shuffling; if omitted, no seeding will be done | |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[Array[String]+]+ | batches | Always | Array of  `nbatches` arrays of the cases listed in  `cases` |
