# In Silico Tasks

These tasks are for use within the bamsurgeon orchestration WDL.

## MutationBED Arrays

In order to introduce genomic variants, bamsurgeon requires the desired mutations to be in a specific BED format. For streamlined use, the bamsurgeon WDL requires a singular input (an parameter called `mutation_bed`) of all desired mutations in a JSON array. The JSON array will be automatically converted to the proper BED format by the bamsurgeon WDL.

All desired mutations in the JSON array must be either all SNVs, SVs, or indels. The JSON array **cannot include a mix of mutation types** (e.g. SNVs with indels, SVs and SNVs, etc). If a mix of these mutation types is desired within a single BAM, run the bamsurgeon WDL to introduce one type of mutation. Then, use the resulting mutated BAM as input to introduce another mutation type.

The information required in the JSON array differs for each mutation type. The following describes what should be included in each type of JSON array:

### SNV JSON Arrays

The inputs requried for SNV mutations includes:

- “chr”: chromosome containing mutation
- “start”: start position of mutation
- “end”: end position of mutation
- “vaf”: variant allele frequency
- “base”: nucleotide to mutate to (optional)

### Indel JSON Arrays

The inputs required for indels includes:

- “chr”: chromosome containing mutation
- “start”: start position of mutation
- “end”: end position of mutation
- “vaf”: variant allele frequency
- "mut_type": either "INS" for insertions or "DEL" for deletions
- "insert_seq": insertion sequence for "INS" mutations (optional)

**Note: When running the workflow to generate indels, it is most reliable to run simple insertions using the "indel" inputs and to run deletions using the "sv" inputs.**

### SV JSON Arrays

SV mutation runs can consist of inversions, deletions, duplications, translocations, and insertions. A mutation will not be made if the largest contig obtained from local assembly of the specified regions is less than three times the maximum library size (max_lib_size). Though small insertions and deletions can be done in sv mutation runs, it is recommended to do an indel run instead.

The required inputs for SV mutations includes:

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

There are additional required inputs for SV translocations:

- “vaf”: variant allele frequency
- “trans_chr”: chromosome for translocation
- “trans_start”: start of translocation
- “trans_end”: end of translocation
- “trans_on_chr”: strand (+ or -) to translocate to
- “trans_from_chr”: strand (+ or -) to translocate from

There are additional required inputs for SV insertions:

- “insert_seq”: insertion sequence
- “vaf”: variant allele frequency
- "target_site_len": length of TSD; for simulating target site duplications (TSDs)
- "cut_site_motif": a sequence of bases with syntax NNN^NNN, where the bases after the caret is the motif; for simulating endocnuclease motifs

For insertions, the sequence to insert can be specified using: a string of the sequence itself, a FASTA file containing the sequence to insert, or an RND library of potential insertions (using the insert_library input). You can simulate target site duplications by including an integer value that specifies the length. Endonuclease motifs can also be simulated by defining an insertion entry with the syntax NNN^NNN, where NNN denotes a sequence of any length, and the cut site motif is denoted by the caret sequence.

There are additional required inputs for SV deletions and inversions:

- “vaf”: variant allele frequency

**Note: When running the workflow to generate indels, it is most reliable to run simple insertions using the "indel" inputs and to run deletions using the "sv" inputs.**

There are additional required inputs for SV duplications:

- “vaf”: variant allele frequency
- “num_dupes”: number of times the sequence should be duplicated

## In Silico Task Parameters

### CramToBamTask Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_cram | Yes | CRAM to convert to BAM | |
| File | input_crai | Yes | Index file for index CRAM | |
| File | ref_fasta | Yes | Reference FASTA file | |
| File | ref_fai | Yes | Index file for reference FASTA | |
| String | output_basename | No | Basename for output BAM files | Defaults to basename of input CRAM |
| String | docker_image | Yes | Docker image with samtools | |
| Int | addldisk | No | Additional disk space to allocate to VM; is added to the disk space needed for the input files | 100 |
| Int | mem_size | No | Amount of memory to allocate to VM | 16 |
| Int | preemptible | No | Number of re-runs allowed for VM | 1 |

### CramToBamTask Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | output_bam | Always | Converted BAM file |
| File | output_bai | Always | Index for converted BAM file |

### RunBamsurgeonTask Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[MutationBED] | mutation_bed | Yes | Array of target mutation information in JSON format | |
| File | input_bam | Yes | BAM file to mutate | |
| File | input_bai | Yes | BAM index for BAM to mutate | |
| File | ref_fasta | Yes | Genome reference file | |
| File | ref_amb | Yes | Genome reference file | |
| File | ref_ann | Yes | Genome reference file | |
| File | ref_bwt | Yes | Genome reference file | |
| File | ref_fai | Yes | Genome reference file | |
| File | ref_pac | Yes | Genome reference file | |
| File | ref_sa | Yes | Genome reference file | |
| File | cnv_file | No | Tabix-indexed list of genome-wide absolute copy number values | |
| File | avoid_reads | No | File of read names to avoid (mutations will be skipped if overlap) | |
| File | donor_bam | No | Bam file for donor reads if using "BIGDUP" mutations | |
| File | donor_bai | No | Bam index file for donor reads if using "BIGDUP" mutations | |
| File | insert_library | No | FASTA file containing library of possible insertions; use "INS RND" instead of "INS filename" to pick one | |
| Int | num_snvs | No | Maximum number of mutations to try | 0 |
| Int | haplo_size | No | Haplotype size | 0 |
| Int | min_depth | No | Minimum read depth to make mutation | 10 |
| Int | max_depth | No | Maximum read depth to make mutation | 2000 |
| Int | min_mut_reads | No | Minimum number of mutated reads to output per site | |
| Int | max_open_files | No | Maximum number of open files during merge | 1000 |
| Int | max_lib_size | No | Maximum fragment length of sequence library | 600 |
| Int | kmer_size | No | Kmer size for assembly | 31 |
| Int | min_contig_gen | No | Minimum length for contig generation, also used to pad assembly | 4000 |
| Int | mean_insert_size | No | Mean insert size | 300 |
| Int | insert_size_stdev | No | Insert size standard deviation | 70 |
| Float | cover_diff | No | Allowable difference in input and output coverage | 0.9 |
| Float | snv_frac | No | Maximum allowable linked SNP MAF (for avoiding haplotypes) | 1 |
| Float | mut_frac | No | Allelic fraction at which to make SNVs | 0.5 |
| Float | sv_frac | No | Allele fraction of variant | 1.0 |
| Float | sim_err | No | Error rate for wgsim-generated reads | |
| Boolean | ignore_snps | No | Make mutations even if there are non-reference alleles sharing the relevant reads | false |
| Boolean | ignore_ref | No | Make mutations even if the mutation is back to the reference allele | false |
| Boolean | ignore_sanity_check | No | Ignore the sanity check enforcing input read count to equal output read count in realignment | false |
| Boolean | single_ended | No | Input bam is single-ended (default is paired-end) | false |
| Boolean | tag_reads | No | Add BS tag to altered reads | false |
| Boolean | ignore_pileup | No | Do not check pileup depth in mutation regions | false |
| Boolean | det_base_changes | No | Deterministic base changes: make transitions only | false |
| Boolean | keep_secondary_reads | No | Keep secondary reads in final BAM | false |
| String | seed | no | Seed for random number generation | |
| String | mutation_type | Yes | Type of mutation to introduce to input BAM file; either "snv", "sv", or "indel" | |
| String | output_basename | Yes | Base name for output files (should include sample name and mutation information) | |
| String | docker_image | Yes | Docker image for running bamsurgeon | |
| Int | addldisk | No | Additional disk space to allocate to VM; is added to the disk space needed for the input files | 200 |
| Int | mem_size | No | Amount of memory to allocate to VM | 8 |
| Int | preemptible | No | Number of re-runs allowed for VM | 1 |

### RunBamsurgeonTask Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | target_bed | Always | Target ROI BED file with mutations to run |
| File | bamsurgeon_script | Always | Bash script that will run bamsurgeon to get desired mutations |
| File | mut_bam | Always | Resulting mutated BAM file from bamsurgeon containing input SNV, indel, or SV mutation(s) |
| File | mut_vcf | Always | VCF containing mutation(s) from bamsurgeon |

### IndexBAMTask Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_bam | Yes | BAM file for sorting and indexing | |
| String | output_basename | No | Basename for output index and sorted BAM | Defaults to the basename of in the input BAM file |
| String | docker_image | Yes | Docker image containing samtools | |
| Int | addldisk | No | Additional disk space to allocate to VM; is added to the disk space needed for the input files | 25 |
| Int | mem_size | No | Amount of memory to allocate to VM | 8 |
| Int | preemptible | No | Number of re-runs allowed for VM | 1 |

### IndexBAMTask Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | output_bam | Always | Sorted BAM file |
| File | output_bai | Always | Index for sorted BAM file |
