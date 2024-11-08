# BGE-PRS Orchestration Workflow

This workflow combines GLIMPSE and PRS WDLs to determine risk scores and percentiles of individuals for a certain disease and model(s).

## GLIMPSE Input Parameters

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

## PRS Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[File] | condition_zip_files | Yes | Array of TAR files, where each TAR contains the variant weights files, score weight file, population loadings file, population meansd file, and population pc file for each condition/disease | |
| File | pruning_sites_for_pca | Yes | Pruning sites used in PCA to generate the ancestry adjustment model | |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |
| String | interaction_docker_image | No | Docker image for use of Python | us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e |
| String | plink_docker_image | No | Docker image for Plink 2 | us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124 |
| String | flash_pca_docker_image | No | Docker image for running flask PCA | us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f |

## Human Genome Reference Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | ref_fasta | Yes | HG38 reference FASTA file | |
| File | ref_fai | Yes | HG38 reference FASTA index file | |
| File | reference_dict | Yes | HG38 Reference dict file | |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | glimpse_vcf | Always | Output imputed VCF from running GLIMPSE |
| File | glimpse_vcf_index | Always | Output imputer VCF index from running GLIMPSE |
| File | glimpse_qc_metrics | When collect_glimpse_qc = "true" | QC metrics from running GLIMPSE |
| Array[File?] | glimpse_phase_monitoring | | |
| File? | glimpse_ligate_monitoring | | |
| File | prs_mix_adjusted_score | Always | Weighted PRS scores adjusted using population models |
| File | pc_projection | Always | PCA projection array |
| File | pc_plot | Always | PCA plot from projection array |
