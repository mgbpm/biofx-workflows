# BGE-PRS Orchestration Workflow

This workflow combines GLIMPSE and PRS WDLs to determine risk scores and percentiles of individuals for a certain disease and model(s).

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | glimpse_reference_chunks | Yes | | |
| Array[File] | input_crams | Yes | Sample input CRAMs to run through GLIMPSE imputation | |
| Array[File] | input_crai | Yes | Sample input CRAM indices to run through GLIMPSE imputation; must be in same order as input_crams | |
| Array[String] | sample_ids | Yes | Sample ids for samples in input_crams | |
| Boolean | impute_reference_only_variants | No | | false |
| Boolean | call_indels | No | | false |
| Int | n_burnin | No | | |
| Int | n_main | No | | |
| Int | effective_population_size | No | | |
| Boolean | collect_glimpse_qc | No | Whether or not to collect qc metrics when running GLIMPSE imputation | true |
| String | glimpse_docker_image | No | Docker image for running GLIMPSE imputation | "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56" |
| String | glimpse_extract_docker_image | No | Docker image to get the number of sites in a GLIMPSE reference chunnk | "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a" |
| File | glimpse_monitoring_script | No | Script to monitor GLIMPSE imputation | |
| Array[File] | condition_files | Yes | Array of TAR files, where each TAR contains the variant weights files, score weight file, population loadings file, population meansd file, and population pc file for each condition/disease | |
| File | scoring_sites | Yes | Sites used in scoring to generate the ancestry adjustment model | |
| File | pruning_sites_for_pca | Yes | Pruning sites used in PCA to generate the ancestry adjustment model | |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |
| File | ref_fasta | Yes | HG38 reference FASTA file | |
| File | ref_fai | Yes | HG38 reference FASTA index file | |
| File | reference_dict | Yes | HG38 Reference dict file | |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | glimpse_vcf | Always | Output imputed VCF from running GLIMPSE |
| File | glimpse_vcf_index | Always | Output imputer VCF index from running GLIMPSE |
| File | glimpse_qc_metrics | When collect_glimpse_qc = "true" | QC metrics from running GLIMPSE |
| Array[File?] | glimpse_phase_monitoring | | |
| File? | glimpse_ligate_monitoring | | |
