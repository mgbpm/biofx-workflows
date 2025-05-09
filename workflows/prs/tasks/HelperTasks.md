## GetBaseMemory Task

Computes the memory (in gigabytes) required by plink to process vcf,
according to the recommendations given in
<https://www.cog-genomics.org/plink/2.0/other#memory>

### Input Parameters

| Type | Name      | Req'd | Description                    | Default Value |
| :--- | :---      | :---  | :---                           | :---          |
| File | vcf       |       | VCF that will be used in PLINK |               |
| Int  | nvariants |       | Number of variants in VCF      |               |

**NB** Exactly one of <code>vcf</code> and <code>nvariants</code> must
be specified.

### Output Parameters

| Type | Name       | When   | Description                     |
| :--- | :---       | :---   | :---                            |
| Int  | gigabytes  | Always | Memory to allocate to PLINK     |
| Int  | nvariants_ | Always | Number of variants in input VCF |

## RenameChromosomesInTsv Task

Renames the human chromosomes mentioned in the first column of
<code>tsv</code> to 1, 2, ..., 22, X, Y, MT.

### Input Parameters

| Type    | Name         | Req'd | Description | Default Value |
| :---    | :---         | :---  | :---        | :---          |
| File    | tsv          | Yes   | TSV with chromosomes to rename | |
| Boolean | skipheader   | Yes   | Whether or not to skip the header of the TSV | |
| File    | lookup       | No    | Chromosome name mapping | "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv" |

### Output Parameters

| Type | Name    | When   | Description                  |
| :--- | :---    | :---   | :---                         |
| File | renamed | Always | TSV with renamed chromosomes |

## RenameChromosomesInVcf Task

Renames the human chromosomes mentioned in the <code>CHROM</code>
column of <code>vcf</code> to 1, 2, ..., 22, X, Y, MT.

### Input Parameters

| Type   | Name         | Req'd | Description | Default Value |
| :---   | :---         | :---  | :---        | :---          |
| File   | vcf          | Yes   | VCF with chromosomes to rename | |
| File   | rename       | No    | Chromosome name mapping | "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv" |

### Output Parameters

| Type | Name    | When   | Description                  |
| :--- | :---    | :---   | :---                         |
| File | renamed | Always | VCF with renamed chromosomes |

## SubsetVcf Task

Subset <code>inputvcf</code> to the regions specified in the
<code>regions</code> file.

### Input Parameters

| Type    | Name      | Req'd | Description                           | Default Value |
| :---    | :---      | :---  | :---                                  | :---          |
| File    | inputvcf  | Yes   | VCF for subsetting                    |               |
| File    | regions   | Yes   | Regions for subsetting the input VCF  |               |
| String  | label     | No    | Basename for the subset VCF           | "data"        |
| Boolean | nocleanup | No    | Include clean-up steps in the process | false         |

### Output Parameters

| Type | Name     | When   | Description              |
| :--- | :---     | :---   | :---                     |
| File | result   | Always | Subsetted VCF            |
| Int  | nregions | Always | Number of regions in VCF |

## Union Task

Generate the set union of the lines in all the files in
<code>lists</code>.  The task simply returns the output of

    sort --unique FILE_0 FILE_1 ... FILE_n

...where <code>FILE_0</code>, <code>FILE_1</code>, ...,
<code>FILE_n</code> are the files in <code>lists</code>

### Input Parameters

| Type         | Name    | Req'd | Description                    | Default Value |
| :---         | :---    | :---  | :---                           | :---          |
| Array[File]+ | lists   | Yes   | (Nonempty) array of text files |               |
| Int          | storage | Yes   | Required storage space, in GB  |               |

### Output Parameters

| Type | Name   | When   | Description                                           |
| :--- | :---   | :---   | :---                                                  |
| File | result | Always | File holding the union of lines in <code>lists</code> |

## Intersection Task

Generate the set intersection of the lines in all the files in
<code>lists</code>.

### Input Parameters

| Type         | Name    | Req'd | Description                    | Default Value |
| :---         | :---    | :---  | :---                           | :---          |
| Array[File]+ | lists   | Yes   | (Nonempty) array of text files |               |
| Int          | storage | Yes   | Required storage space, in GB  |               |

### Output Parameters

| Type | Name   | When   | Description                                               |
| :--- | :---   | :---   | :---                                                      |
| File | result | Always | File holding the intersection of lines <code>lists</code> |

## ListShards Task

List all the VCF shard files in the directory tree rooted at
<code>source</code>.  (Here, "VCF shard files" means files with
extension ".vcf.gz" whose parent directory's basename is "shards".)

### Input Parameters

| Type   | Name         | Req'd | Description                            | Default Value |
| :---   | :---         | :---  | :---                                   | :---          |
| String | source       | Yes   | URL to a GCP bucket                    |               |
| String | workspace    | Yes   | name of workspace where Terra job runs |               |
| String | docker_image | Yes   | name of docker image for the job       |               |

### Output Parameters

| Type | Name     | When   | Description                                                           |
| :--- | :---     | :---   | :---                                                                  |
| File | relpaths | Always | File the paths of all VCF shard paths relative to <code>source</code> |

## MakeBatches Task

Distribute the items in <code>cases</code> across
<code>nbatches</code> batches of roughly equal size.

### Input Parameters

| Type          | Name      | Req'd | Description                                            | Default Value |
| :---          | :---      | :---  | :---                                                   | :---          |
| File          | cases     | Yes   | File containing cases to batch, one per line           |               |
| Int           | nbatches  | Yes   | Desired number of batches                              |               |
| Array[String] | exclude   | No    | (Possibly empty) array of exclusions                   | []            |
| Boolean       | noshuffle | No    | If true, do not shuffle the cases before batching them | false         |
| Int?          | seed      | No    | Seed to use for the PRNG before shuffling; if omitted, no seeding will be done | |

### Output Parameters

| Type                   | Name    | When   | Description                                                                     |
| :---                   | :---    | :---   | :---                                                                            |
| Array[Array[String]+]+ | batches | Always | Array of <code>nbatches</code> arrays of the cases listed in <code>cases</code> |
