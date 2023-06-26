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
| String | assay\_type | yes | Either "amplicon" or "capture" (for Hybrid capture) | |
| String | docker\_image | yes | name of insiM docker image to use | |
| File | genome\_ref\_sa | yes | genome reference file found in workspace | |
| File | genome\_ref\_amb | yes | genome reference file found in workspace | |
| File | genome\_ref\_ann | yes | genome reference file found in workspace | |
| File | genome\_ref\_bwt | yes | genome reference file found in workspace | |
| File | genome\_ref\_fasta | yes | genome reference file found in workspace | |
| File | genome\_ref\_pac | yes | genome reference file found in workspace | |
| File | input\_bai | yes | index file corresponding to bam file needing mutation | |
| File | input\_bam | yes | bam file needing mutation; generally a NIST sample located in the LMM reference data | |
| String | mutation\_type | yes | desired type of mutation; options are "snv", "ins", "del", "dup", or "mix" | |
| String | output\_basename | yes | desired name for mutated bam; ex: "mut\_test" will make mut\_test.bam | |
| Array | targets | yes | Json array of target regions to mutate | |
| String | amp\_read\_length | no | read length (only necessary when assay\_type = "amplicon") | |
| Array | ampbed | no | Json array including amplicon mutations; should look like the targets array (only necessary when assay_type = "amplicon") | |
| String | inslen | no | size/length of insert | 10 |
| String | insseq | no | sequence of insert | a random sequence of length "inslen" |
| String | vaf | no | variant allele frequency (only necessary for mutations not of type "mix") | |


<u>Notes on Parameters</u>
For 'mix' mutations:
- indicate the types of individual mutations within the targets array 'type'
- indicate the size/length of insert for individual mutations within the targets array 'ins\_len'
- indicate the sequence of insert for individual mutations within the targets array 'ins\_seq'

The targets array should include the following information for each mutation to be included:
- “chr”: chromosome containing mutation
- “start”: start position of mutation
- “end”: end position of mutation
- “type”: type of mutation; either 'snv', 'ins', 'del', or 'dup'
- “vaf”: variant allele frequency
- “ins\_seq”: insert sequence (optional)
- “ins\_len”: length of insert (optional)

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | fastq1 | Always | first mutated FASTQ based off the input bam; includes mutation from targets array |
| File | fastq2 | Always | second mutated FASTQ based off the input bam; includes mutation from targets array |
| File | vcf | Always | VCF file corresponding to mutated FASTQs |
| File | mut\_sam | Always | SAM file produced from aligning mutated FASTQs to reference genome using bwa-mem |
| File | mut\_sorted\_bam | Always | BAM file produced by sorting mut\_sam |
| File | mut\_dedup\_bam | Always | BAM file produced by marking/removing duplications in mut\_sorted\_bam |
| File | metrics\_file | Always | text file with metrics from deduplication of mut\_sorted\_bam |
| File | mut\_dedup\_bai | Always | index file for mut\_dedup\_bam |