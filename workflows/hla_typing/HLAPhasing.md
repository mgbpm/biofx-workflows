# HLA Phasing Workflow

This workflow will run phasing on the output of the T1K Genotyping WDL, specifically the allele TSV file. All source code and config files for running phasing can be found in the biofx-prs repository under the hla_imputation directory.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_alleles | Yes | TSV of HLA alleles, which is the output of the allele TSV from T1K; each allele must be on a separate line | |
| String | prs_docker_image | No | Docker image from biofx-prs repository, which includes the source code and configs for phasing | |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | scored_haplotypes | Always | Text file with haplotype codes, weights, and final T1D_MAPS score for the sample |
