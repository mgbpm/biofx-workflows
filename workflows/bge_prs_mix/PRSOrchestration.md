# BGE-PRS Orchestration Workflow

This workflow combines GLIMPSE and PRS WDLs to determine risk scores and percentiles of individuals for a certain disease and model(s).

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | ref_fasta | Yes | HG38 reference FASTA file | |
| File | ref_fai | Yes | HG38 reference FASTA index file | |
| File | reference_dict | Yes | HG38 Reference dict file | |
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
| String | glimpse_extrat_docker_image | No | Docker image to get the number of sites in a GLIMPSE reference chunnk | "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a" |
| File | glimpse_monitoring_script | No | Script to monitor GLIMPSE imputation | |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | glimpse_vcf | Always | Output imputed VCF from running GLIMPSE |
| File | glimpse_vcf_index | Always | Output imputer VCF index from running GLIMPSE |
| File | glimpse_qc_metrics | When collect_glimpse_qc = "true" | QC metrics from running GLIMPSE |
| Array[File?] | glimpse_phase_monitoring | | |
| File? | glimpse_ligate_monitoring | | |
