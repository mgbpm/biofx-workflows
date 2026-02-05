# Running a bamsurgeon workflow

## Background
<u>Usage</u>
bamsurgeon is a tool that will introduce genomic variants (such as SNV, SV, and indel mutations) into BAM, CRAM, and SAM files. The resulting mutated BAM files can be used to test variant callers.

<u>About the Workflow</u>
Bamsurgeon was originally developed by Adam Ewing. This workflow adapts the original Dockerfile and python scripts to work with a WDL on Terra. All scripts that are used in the tasks of this workflow can be found on the MGBPM-IT Azure DevOps Repo for bamsurgeon.

bamsurgeon has three main scripts, each corresponding to three categories of mutations that can be introduced to BAM files: addsnv.py, addindel.py, and addsv.py. Each of the different mutations has various optional inputs to run the workflow, most of which have defaults listed in the tables under the “Input Variables” section below.

## Input Parameters for all Mutation Types
Note that BWA makes assumptions that all genome reference files are conventionally named. The reference files found under the hg38 reference data in Terra suit this assumption. In addition, these reference files are not explicitly used in the WDL (besides the FASTA), but are necessary within the container for BWA to run.

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[MutationBED] | mutation_bed_input | Yes | Array of target mutation information in JSON format | | 
| String | output_files_base | Yes | Base name for output files (should include sample name and mutation information) | |
| String | mutation_type | Yes | Type of mutation to introduce to input BAM file; either "snv", "sv", or "indel" | |
| File | input_bam_file | Yes | BAM file to mutate | |
| File | input_bai_file | Yes | BAM index for BAM to mutate | |
| File | ref_fasta | Yes | Genome reference file found in workspace | |
| File | ref_amb | Yes | Genome reference file found in workspace | |
| File | ref_ann | Yes | Genome reference file found in workspace | |
| File | ref_bwt | Yes | Genome reference file found in workspace | |
| File | ref_fai | Yes | Genome reference file found in workspace | |
| File | ref_pac | Yes | Genome reference file found in workspace | |
| File | ref_sa | Yes | Genome reference file found in workspace | |
| File | ref_dict | No | Genome reference file found in workspace; required for PGx/Risk | |
| File | dbsnp_vcf | No | dbSNP VCF file that matches the genome reference; required for PGx/Risk | |
| File | dbsnp_vcf_index | No | dbSNP VCF file index; required for PGx/Risk | |
| String | bamsurgeon_docker_image | No | Name of bamsurgeon docker image | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bamsurgeon:20240229" |
| String | samtools_docker_image | No | Name of samtools docker image | "biocontainers/samtools:v1.9-4-deb_cv1" |
| String | igvreport_docker_image | No | Name of IGV report docker image | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511" |
| Boolean | run_pgx | No | Whether or not to run the PGx pipeline with the mutated BAM | true |
| Boolean | run_risk | No | Whether or not to run the Risk pipeline with the mutated BAM | true |
| String | pgx_test_code | No | Test code that defines which pharmacogenomics report to generate | "lmPGX-pnlD_L" |
| String | pgx_docker_image | No | The name of the Docker image to generate the pharmacogenomics report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20241007" |
| File | pgx_workflow_fileset | No | Tar file containing the pharmacogenomics reference data to generate the report | "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_20241004.tar" |
| File | pgx_roi_bed | No | BED file that defines the genomic regions to include in the pharmacogenomics analysis | "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_genotyping.bed" |
| String | risk_alleles_test_code | No | Test code that defines which risk alleles report to generate | "lmRISK-pnlB_L" |
| String | risk_alleles_docker_image | No | The name of the Docker image to generate the risk alleles report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20240129" |
| File | risk_alleles_workflow_fileset | No | Tar file containing the risk alleles reference data to generate the report | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar" |
| File | risk_alleles_roi_bed | No | BED file that defines the genomic regions to include in the risk alleles analysis | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed" |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | No | The name of the current workspace (for secret retrieval); require for running PGx/Risk | |

## Input Parameters for SNV-Specific Mutations
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Float | snv_frac | No | Maximum allowable linked SNP MAF (for avoiding haplotypes) | 1 |
| Float | mut_frac | No | Allelic fraction at which to make SNVs | 0.5 |
| Int | num_snvs | No | Maximum number of mutations to try | 0 |
| File | cnv_file | No | Tabix-indexed list of genome-wide absolute copy number values | |
| Float | cover_diff | No | Allowable difference in input and output coverage | 0.9 |
| Int | haplo_size | No | Haplotype size | 0 |
| Int | min_depth | No | Minimum read depth to make mutation | 10 |
| Int | max_depth | No | Maximum read depth to make mutation | 2000 |
| Int | min_mut_reads | No | Minimum number of mutated reads to output per site | |
| File | avoid_reads | No | File of read names to avoid (mutations will be skipped if overlap) | |
| Boolean | ignore_snps | No | Make mutations even if there are non-reference alleles sharing the relevant reads | false |
| Boolean | ignore_ref | No | Make mutations even if the mutation is back to the reference allele | false |
| Boolean | ignore_sanity_check | No | Ignore the sanity check enforcing input read count to equal output read count in realignment | false |
| Boolean | single_ended | No | Input bam is single-ended (default is paired-end) | false |
| Int | max_open_files | No | Maximum number of open files during merge | 1000 |
| Boolean | tag_reads | No | Add BS tag to altered reads | false |
| Boolean | ignore_pileup | No | Do not check pileup depth in mutation regions | false |
| String | seed | no | Seed for random number generation | |

The information to be included in the snv_bed_input array includes:
- “chr”: chromosome containing mutation
- “start”: start position of mutation
- “end”: end position of mutation
- “vaf”: variant allele frequency
- “base”: nucleotide to mutate to (optional)


## Input Parameters for Indel-Specific Mutations
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Float | snv_frac | No | Maximum allowable linked SNP MAF (for avoiding haplotypes) | 1 |
| Float | mut_frac | No | Allelic fraction at which to make SNVs | 0.5 |
| Int | num_snvs | No | Maximum number of mutations to try | 0 |
| File | cnv_file | No | Tabix-indexed list of genome-wide absolute copy number values | |
| Float | cover_diff | No | Allowable difference in input and output coverage | 0.9 |
| Int | min_depth | No | Minimum read depth to make mutation | 10 |
| Int | max_depth | No | Maximum read depth to make mutation | 2000 |
| Int | min_mut_reads | No | Minimum number of mutated reads to output per site | |
| File | avoid_reads | No | File of read names to avoid (mutations will be skipped if overlap) | |
| Boolean | det_base_changes | No | Deterministic base changes: make transitions only | false |
| Boolean | ignore_sanity_check | No | Ignore the sanity check enforcing input read count to equal output read count in realignment | false |
| Boolean | single_ended | No | Input bam is single-ended (default is paired-end) | false |
| Int | max_open_files | No | Maximum number of open files during merge | 1000 |
| Boolean | tag_reads | No | Add BS tag to altered reads | false |
| Boolean | ignore_pileup | No | Do not check pileup depth in mutation regions | false |
| String | seed | No | See for random number generation | |

The information to be included in the indel_bed_input array includes:
- “chr”: chromosome containing mutation
- “start”: start position of mutation
- “end”: end position of mutation
- “vaf”: variant allele frequency
- "mut_type": either "INS" for insertions or "DEL" for deletions
- "insert_seq": insertion sequence for "INS" mutations (optional)

**Note: When running the workflow to generate indels, it is most reliable to run simple insertions using the "indel" inputs and to run deletions using the "sv" inputs.**


## Input Parameters for SV-Specific Mutations
SV mutation runs can consist of inversions, deletions, duplications, translocations, and insertions. A mutation will not be made if the largest contig obtained from local assembly of the specified regions is less than three times the maximum library size (max_lib_size). Though small insertions and deletions can be done in sv mutation runs, it is recommended to do an indel run instead.

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Int | max_lib_size | No | Maximum fragment length of sequence library | 600 |
| Int | kmer_size | No | Kmer size for assembly | 31 |
| Float | sv_frac | No | Allele fraction of variant | 1.0 |
| Int | min_contig_gen | No | Minimum length for contig generation, also used to pad assembly | 4000 |
| File | cnv_file | No | Tabix-indexed list of genome-wide absolute copy number values | |
| File | donor_bam | No | Bam file for donor reads if using "BIGDUP" mutations | |
| File | donor_bai | No | Bam index file for donor reads if using "BIGDUP" mutations | |
| Int | mean_insert_size | No | Mean insert size | 300 |
| Int | insert_size_stdev | No| Insert size standard deviation | 70 |
| Float | sim_err | No | Error rate for wgsim-generated reads | |
| File | insert_library | No | FASTA file containing library of possible insertions; use "INS RND" instead of "INS filename" to pick one | |
| Boolean | tag_reads | No | Add BS tag to altered reads | false |
| Boolean | keep_secondary_reads | No | Keep secondary reads in final BAM | false |
| String | seed | No | See for random number generation | |

The information in the sv_bed_input array should include:
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

Additional information in sv_input array for translocations:
- “vaf”: variant allele frequency
- “trans_chr”: chromosome for translocation
- “trans_start”: start of translocation
- “trans_end”: end of translocation
- “trans_on_chr”: strand (+ or -) to translocate to
- “trans_from_chr”: strand (+ or -) to translocate from

Additional information in sv_input array for insertions:
- “insert_seq”: insertion sequence
- “vaf”: variant allele frequency
- "target_site_len": length of TSD; for simulating target site duplications (TSDs)
- "cut_site_motif": a sequence of bases with syntax NNN^NNN, where the bases after the caret is the motif; for simulating endocnuclease motifs

For insertions, the sequence to insert can be specified using: a string of the sequence itself, a FASTA file containing the sequence to insert, or an RND library of potential insertions (using the insert_library input). You can simulate target site duplications by including an integer value that specifies the length. Endonuclease motifs can also be simulated by defining an insertion entry with the syntax NNN^NNN, where NNN denotes a sequence of any length, and the cut site motif is denoted by the caret sequence.

Additional information in sv_input array for deletions and inversions:
- “vaf”: variant allele frequency

**Note: When running the workflow to generate indels, it is most reliable to run simple insertions using the "indel" inputs and to run deletions using the "sv" inputs.**

Additional information in sv_input for duplications:
- “vaf”: variant allele frequency
- “num_dupes”: number of times the sequence should be duplicated

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | target_bed | Always | Target ROI BED file with mutations to run |
| File | bamsurgeon_script | Always | Bash script that will run bamsurgeon to get desired mutations |
| File | mutated_bam | Always | Resulting sorted BAM file from bamsurgeon containing input SNV, indel, or SV mutation |
| File | mutated_bai | Always | Index file corresponding to the final sorted BAM file |
| File | mutated_vcf | Always | VCF containing mutations from bamsurgeon |
| File | igv_report | Always | IGV report detailing desired mutations in final BAM |
| File | pgx_details_report | Always | PGx detailed Excel report with Drug recommendations|
| File | pgx_summary_report | Always | PGx summarized Excel report with links to publications |
| File | pgx_genotype_xlsx | If PGx is enabled | Full list of pharmacogenomics genotypes in XLSX format |
| File | pgx_genotype_txt | If PGx is enabled | Full list of pharmacogenomics genotypes in TSV format |
| File | risk_alleles_report | If Risk is enabled | Risk alleles report |
| File | risk_alleles_genotype_xlsx | If Risk is enabled | Full list of risk allele genotypes in XLSX format |
| File | risk_alleles_genotype_txt | If Risk is enabled | Full list of risk allele genotypes in TSV format |
