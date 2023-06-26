# Running a bamsurgeon workflow

## Background
<u>Usage</u>
bamsurgeon is a tool that will introduce genomic variants (such as SNV, SV, and indel mutations) into BAM, CRAM, and SAM files. The resulting mutated BAM files can be used to test variant callers.

<u>About the Workflow</u>
Bamsurgeon was originally developed by Adam Ewing. This workflow adapts the original Dockerfile and python scripts to work with a WDL on Terra. All scripts that are used in the tasks of this workflow can be found on the[MGBPM-IT Azure DevOps Repo for bamsurgeon.

bamsurgeon has three main scripts, each corresponding to three categories of mutations that can be introduced to BAM files: addsnv.py, addindel.py, and addsv.py. Each of the different mutations has various optional inputs to run the workflow, most of which have defaults listed in the tables under the “Input Variables” section below.

## Input Parameters
<u>Mandatory Variables for all Mutation Types</u>
Note that BWA makes assumptions that all genome reference files are conventionally named. The reference files found under the hg38 reference data in Terra suit this assumption. In addition, these reference files are not explicitly used in the WDL (besides the FASTA), but are necessary within the container for BWA to run.

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | mutation\_to\_run | yes | Type of mutation to be run; either "snv", "sv", or "indel" | |
| String | out\_bam\_name | yes | Bam file name for the output | |
| String | docker\_image | yes | Name of docker image to run WDL with | |
| File | bai\_file | yes | index for bam file that will be mutated | |
| File | bam\_file | yes | bam file that will be mutate | |
| File | genome\_ref\_sa | yes | genome reference file found in workspace | |
| File | genome\_ref\_amb | yes | genome reference file found in workspace | |
| File | genome\_ref\_ann | yes | genome reference file found in workspace | |
| File | genome\_ref\_bwt | yes | genome reference file found in workspace | |
| File | genome\_ref\_fasta | yes | genome reference file found in workspace | |
| File | genome\_ref\_pac | yes | genome reference file found in workspace | |


<u>Variables Specific to Running SNV Mutations</u>

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array | snv\_input | yes | Array of target mutation information in JSON format |  |
| Float | snv\_frac | no | Maximum allowable linked SNP MAF (for avoiding haplotypes) | 1 |
| Float | mut\_frac | no | Allelic fraction at which to make SNVs | 0.5 |
| Int | num\_snvs | no | Maximum number of mutations to try | 0 |
| File | cnv\_file | no | Tabix-indexed list of genome-wide absolute copy number values | |
| Float | cover\_diff | no | Allowable difference in input and output coverage | 0.9 |
| Int | haplo\_size | no | Haplotype size | 0 |
| Int | min\_depth | no | Minimum read depth to make mutation | 10 |
| Int | max\_depth | no | Maximum read depth to make mutation | 2000 |
| Int | min\_mut\_reads | no | Minimum number of mutated reads to output per site | |
| File | avoid\_reads | no | File of read names to avoid (mutations will be skipped if overlap) | |
| Boolean | ignore\_snps | no | Make mutations even if there are non-reference alleles sharing the relevant reads | False |
| Boolean | ignore\_ref | no | Make mutations even if the mutation is back to the reference allele | False |
| Boolean | ignore\_sanity\_check | no | Ignore the sanity check enforcing input read count to equal output read count in realignment | False |
| Boolean | single\_ended | no | Input bam is single-ended (default is paired-end) | False |
| Int | max\_open\_files | no | Maximum number of open files during merge | 1000 |
| Boolean | tag\_reads | no | Add BS tag to altered reads | False |
| Boolean | ignore\_pileup | no | Do not check pileup depth in mutation regions | False |
| String | seed | no | See for random number generation | |

The information to be included in the snv_input array includes:
- “chr”: chromosome containing mutation
- “start”: start position of mutation
- “end”: end position of mutation
- “vaf”: variant allele frequency
- “base”: nucleotide to mutate to (optional)


<u>Variables Specific to Running Indel Mutations</u>

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array | indel\_input | yes | Array of target mutation information in JSON format |  |
| Float | snv\_frac | no | Maximum allowable linked SNP MAF (for avoiding haplotypes) | 1 |
| Float | mut\_frac | no | Allelic fraction at which to make SNVs | 0.5 |
| Int | num\_snvs | no | Maximum number of mutations to try | 0 |
| File | cnv\_file | no | Tabix-indexed list of genome-wide absolute copy number values | |
| Float | cover\_diff | no | Allowable difference in input and output coverage | 0.9 |
| Int | min\_depth | no | Minimum read depth to make mutation | 10 |
| Int | max\_depth | no | Maximum read depth to make mutation | 2000 |
| Int | min\_mut\_reads | no | Minimum number of mutated reads to output per site | |
| File | avoid\_reads | no | File of read names to avoid (mutations will be skipped if overlap) | |
| Boolean | deterministic\_base\_changes | no | Deterministic base changes: make transitions only | False |
| Boolean | ignore\_sanity\_check | no | Ignore the sanity check enforcing input read count to equal output read count in realignment | False |
| Boolean | single\_ended | no | Input bam is single-ended (default is paired-end) | False |
| Int | max\_open\_files | no | Maximum number of open files during merge | 1000 |
| Boolean | tag\_reads | no | Add BS tag to altered reads | False |
| Boolean | ignore\_pileup | no | Do not check pileup depth in mutation regions | False |
| String | seed | no | See for random number generation | |

The information to be included in the snv\_input array includes:
- “chr”: chromosome containing mutation
- “start”: start position of mutation
- “end”: end position of mutation
- “vaf”: variant allele frequency
- "mut\_type": either "INS" for insertions or "DEL" for deletions
- "insert\_seq": insertion sequence for "INS" mutations (optional)


<u>Variables Specific to Running SV Mutations</u>
SV mutation runs can consist of inversions, deletions, duplications, translocations, and insertions. A mutation will not be made if the largest contig obtained from local assembly of the specified regions is less than three times the maximum library size (max_lib_size). Though small insertions and deletions can be done in sv mutation runs, it is recommended to do an indel run instead.

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array | sv\_input | yes | Array of target mutation information in JSON format |  |
| Int | max\_lib\_size | no | Maximum fragment length of sequence library | 600 |
| Int | kmer\_size | no | Kmer size for assembly | 31 |
| Float | sv\_frac | no | Allele fraction of variant | 1.0 |
| Int | min\_contig\_gen | no | Minimum length for contig generation, also used to pad assembly | 4000 |
| File | cnv\_file | no | Tabix-indexed list of genome-wide absolute copy number values | |
| File | donor\_bam | no | Bam file for donor reads if using "BIGDUP" mutations | |
| File | donor\_bai | no | Bam index file for donor reads if using "BIGDUP" mutations | |
| Int | mean\_insert\_size | no | Mean insert size | 300 |
| Int | insert\_size\_stdev | no| Insert size standard deviation | 70 |
| Float | sim\_err | no | Error rate for wgsim-generated reads | |
| File | insert\_library | no | FASTA file containing library of possible insertions; use "INS RND" instead of "INS filename" to pick one | |
| Boolean | tag\_reads | no | Add BS tag to altered reads | False |
| Boolean | keep\_secondary\_reads | no | Keep secondary reads in final BAM | False |
| String | seed | no | See for random number generation | |

The information in the sv\_input array should include:
- “chr”: chromosome containing mutation
- “start”: start position of mutation
- “end”: end position of mutation
- “mut_type":
    - "INS" for insertions
	- “DEL” for deletions
	- “BIGDEL” for deletions greater than 5 kbp in length
	- “DUP” for duplications
	- “BIGDUP” for duplications greater than 5 kbp
	- “INV” for inversions
	- “BIGINV” for inversions greater than 5 kbp in length
	- “TRN” for translocation

Additional information in sv\_input array for translocations:
- “vaf”: variant allele frequency
- “trans\_chr”: chromosome for translocation
- “trans\_start”: start of translocation
- “trans\_end”: end of translocation
- “trans\_on\_chr”: strand (+ or -) to translocate to
- “trans\_from\_chr”: strand (+ or -) to translocate from

Additional information in sv\_input array for insertions:
- “insert\_seq”: insertion sequence
- “vaf”: variant allele frequency
- "target\_site\_len": length of TSD; for simulating target site duplications (TSDs)
- "cut\_site\_motif": a sequence of bases with syntax NNN^NNN, where the bases after the caret is the motif; for simulating endocnuclease motifs

For insertions, the sequence to insert can be specified using: a string of the sequence itself, a FASTA file containing the sequence to insert, or an RND library of potential insertions (using the insert_library input). You can simulate target site duplications by including an integer value that specifies the length. Endonuclease motifs can also be simulated by defining an insertion entry with the syntax NNN^NNN, where NNN denotes a sequence of any length, and the cut site motif is denoted by the caret sequence.

Additional information in sv_input array for deletions and inversions:
- “vaf”: variant allele frequency

Additional information in sv\_input for duplications:
- “vaf”: variant allele frequency
- “num\_dupes”: number of times the sequence should be duplicated