# Running an insiM workflow

## Background
<u>Usage</u>
insiM (in silico Mutator) is a tool that will introduce point mutations, insertions, deletions, and duplications of any size into real datasets of amplicon-based or hybrid-capture NGS assay (Patil, 2019). The resulting mutated files can be used to validated bioinformatics pipelines, such as variant calling rare mutations.

Citation: Patil, Sushant A., et al. 'insiM: in silico mutator software for bioinformatics pipeline validation of clinical next-generation sequencing assays.' The Journal of Molecular Diagnostics 21.1 (2019): 19-26.

<u>About the Workflow</u>
insiM was originally developed by Sushant Patil and team to run locally (not in WDL form) and with Python 2.7.6. This insiM workflow adapted the original insiM scripts to use Python 3.11. All scripts that are used in tasks within the workflow can be found on the MGBPM-IT Azure DevOps Repo for insiM.

insiM, by design, takes in a BAM file for mutation and produces paired-end FASTQ files containing desired mutations. This workflow will also align these FASTQ files with reference files using BWA-mem. It will also convert the resulting SAM file to a final BAM (without duplicate reads). An index file for the final BAM will also be produced.

## Input Parameters
Note that BWA makes assumptions that all genome reference files are conventionally named. The reference files found under the hg38 reference data in Terra suit this assumption. In addition, these reference files are not explicitly used in the WDL (besides the FASTA), but are necessary within the container for BWA to run.

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_bam | Yes | BAM file needing mutation; generally a NIST sample located in the LMM reference data | |
| File | input_bai | Yes | Index file corresponding to BAM file needing mutation | |
| File | ref_fasta | Yes | Genome reference file found in workspace | |
| File | ref_fai | Yes | Genome reference file found in workspace | |
| File | ref_amb | Yes | Genome reference file found in workspace | |
| File | ref_ann | Yes | Genome reference file found in workspace | |
| File | ref_bwt | Yes | Genome reference file found in workspace | |
| File | ref_pac | Yes | Genome reference file found in workspace | |
| File | ref_sa | Yes | Genome reference file found in workspace | |
| File | ref_dict | No | Genome reference file found in workspace; required for PGx/Risk | |
| File | dbsnp_vcf | No | dbSNP VCF file that matches the genome reference; required for PGx/Risk | |
| File | dbsnp_vcf_index | No | dbSNP VCF file index; required for PGx/Risk | |
| Array | targets | Yes | JSON array of target regions to mutate; includes chr, start, end, type, vaf, insert sequence, and insert sequence length | |
| Array | ampbed | No | JSON array including amplicon mutations; should look like the targets array (only necessary when assay_type = "amplicon") | |
| String | assay_type | Yes | Either "amplicon" or "capture" (for Hybrid capture) | |
| String | mutation_type | Yes | Desired type of mutation; options are "snv", "ins", "del", "dup", or "mix" | |
| String | vaf | No | Variant allele frequency (only necessary for mutations not of type "mix") | |
| String | inslen | No | Size/length of insert | No default is given in the WDL, but insiM by default uses a length of 10 |
| String | insseq | No | sequence of insert | No default is given in the WDL, but insiM by default uses a random sequence of length "inslen" |
| String | amp_read_length | No | Read length (only necessary when assay_type = "amplicon") | |
| String | sample_name | Yes | Name of sample | |
| String | mutation_name | Yes | Name for intended mutation | |
| String | insim_docker_image | No | Name of insiM docker image to use | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/insim:testing" |
| String | picard_docker_image | No | Name of picard docker image to use | "biocontainers/picard:v1.139_cv3" |
| String | bwa_docker_image | No | Name of BWA docker image to use | "biocontainers/bwa:v0.7.17_cv1" |
| String | samtools_docker_image | No | Name of samtools docker image to use | "biocontainers/samtools:v1.9-4-deb_cv1" |
| String | igvreport_docker_image | No | Name of IGV report docker image | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511" |
| Boolean | run_pgx | No | Whether or not to run the PGx pipeline with the mutated BAM | true |
| Boolean run_risk | No | Whether or not to run the Risk pipeline with the mutated BAM | true |
| String | pgx_test_code | No | Test code that defines which pharmacogenomics report to generate | "lmPGX-pnlC_L" |
| String | pgx_docker_image | No | The name of the Docker image to generate the pharmacogenomics report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20240129" |
| File | pgx_workflow_fileset | No | Tar file containing the pharmacogenomics reference data to generate the report | "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_files-20220118.tar" |
| File | pgx_roi_bed | No | BED file that defines the genomic regions to include in the pharmacogenomics analysis | "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_genotyping-chr-20220118.bed" |
| String | risk_alleles_test_code | No | Test code that defines which risk alleles report to generate | "lmRISK-pnlB_L" |
| String | risk_alleles_docker_image | No | The name of the Docker image to generate the risk alleles report | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20240129" |
| File | risk_alleles_workflow_fileset | No | Tar file containing the risk alleles reference data to generate the report | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar" |
| File | risk_alleles_roi_bed | No | BED file that defines the genomic regions to include in the risk alleles analysis | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed" |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | No | The name of the current workspace (for secret retrieval); require for running PGx/Risk | |


<u>Notes on Parameters</u>
For 'mix' mutations:
- indicate the types of individual mutations within the targets array 'type'
- indicate the size/length of insert for individual mutations within the targets array 'ins_len'
- indicate the sequence of insert for individual mutations within the targets array 'ins_seq'

The targets array should include the following information for each mutation to be included:
- “chr”: chromosome containing mutation
- “start”: start position of mutation
- “end”: end position of mutation
- “type”: type of mutation; either 'snv', 'ins', 'del', or 'dup'
- “vaf”: variant allele frequency
- “ins_seq”: insert sequence (optional)
- “ins_len”: length of insert (optional)

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | output_fastq1 | Always | First mutated FASTQ based off the input BAM; includes mutation from targets array |
| File | output_fastq2 | Always | Second mutated FASTQ based off the input BAM; includes mutation from targets array |
| File | output_vcf | Always | VCF file corresponding to mutated FASTQs |
| File | output_mutated_sam | Always | SAM file produced from aligning mutated FASTQs to reference genome using bwa-mem |
| File | output_dedup_bam | Always | BAM file produced by marking/removing duplications in sorted mutated BAM |
| File | output_dedup_metrics | Always | TXT file with metrics from deduplication |
| File | output_dedup_bai | Always | Index file for output_dedup_bam |