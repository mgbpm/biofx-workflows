## GetBaseMemory Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | vcf | No | VCF that will be used in PLINK | |
| Int | nvariants | No | Number of variants in VCF | |
| String | docker_image | No | Python Docker image | "python:3.11" |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| Int | gigabytes | Always | Memory to allocate to PLINK |
| Int | nvariants_ | Always | Number of variants in input VCF |

## RenameChromosomesInTsv Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | tsv | Yes | TSV with chromosomes to rename | |
| Boolean | skipheader | Yes | Whether or not to skip the header of the TSV | |
| File | lookup | No | Chromosome renaming file for conventions | "gs://fc-secure-9ea53c3d-d71a-4f59-92c3-63c75c622a88/reference/etc/rename_chromosomes.tsv" |
| String | docker_image | No | Python Docker image | "python:3.11" |
| Int | disk_size | No | Disk size to allocate in GB | Size of input VCF (x2) plus 20 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | renamed | Always | TSV with renamed chromosomes |

## RenameChromosomesInVcf Task

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | vcf | Yes | VCF with chromosomes to rename | |
| File | rename | No | Chromosome renaming file for conventions | "gs://fc-secure-9ea53c3d-d71a-4f59-92c3-63c75c622a88/reference/etc/rename_chromosomes.tsv" |
| String | docker_image | No | Python Docker image | "biocontainers/bcftools:v1.9-1-deb_cv1" |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | renamed | Always | VCF with renamed chromosomes |

## SubsetVcf

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | inputvcf | Yes | VCF for subsetting | |
| File | regions | Yes | Regions for subsetting the input VCF | |
| String | label | No | Label for the data | "data" |
| Boolean | nocleanup | No | Implement clean-up process | false |
| String | docker_image | No | Docker image for bcftools | "biocontainers/bcftools:v1.9-1-deb_cv1" |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | result | Always | Subsetted VCF |
| Int | nregions | Always | Number of regions in VCF |