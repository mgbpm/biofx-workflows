# ReblockGVCFandPush Workflow
workflow assumes input will be a gvcf file located in a local GCP bucket
workflow reblocks input gvcf files producting reblocked gvcf files as output
Reblocking perfromed using the Broad's RebockGVCFv2.2.1 WDL
Found at: "https://raw.githubusercontent.com/broadinstitute/warp/f0e6d797fef941c2cfea260a9dd0adcb8effe961/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl"

## Input Parameters
input_gvcf: gvcf file to reblock
input_gvcf_index: index file of gvcf to reblock
rbgvcf_staging_bucket: gcp bucket location to transfer reblocked gvcf file + index to

## Outputs
Reblock.output_vcf: Reblocked gvcf file
Reblock.output_vcf_index: Index file for Reblocked gvcf file