# PRS Mix Orchestration Workflow

This workflow combines GLIMPSE and PRS WDLs to determine risk scores and percentiles of individuals for a certain disease and model(s). The workflow is designed for single sample input. The model manifests used for each desired condition in the assay are to be generated before running this workflow using the MakeAdjustmentModelWorkflow WDL. It is also recommended that the inputs for running the model -- and thus those consequently used in this workflow -- are cleaned using the PreparePrsMixInputsWorkflow WDL to ensure no errors occur when using flashpca.

This workflow is also designed to be used with Sample Tracker.

## Fetch File Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | data_location | Yes | Source location of CRAM and CRAI | |
| String | sample_id | Yes | Sample ID for CRAM and CRAI | |
| String | subject_id | Yes | Subject ID to match the input sample ID | |
| String | orchutils_docker_image | No | Docker image for orchestration tasks, such as fetching CRAM and CRAI | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20250203" |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| Int | fetch_disk_size | No | Size of disk to allocation in GB for fetching CRAM and CRAI | 75 |

## GLIMPSE Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | glimpse_reference_chunks | Yes | List of file paths to files that contain reference chunks for GLIMPSE | |
| File | ref_fasta | Yes | HG38 reference FASTA file | |
| File | ref_fai | Yes | HG38 reference FASTA index file | |
| File | reference_dict | Yes | HG38 Reference dict file | |
| File | gnomadAF_ref_vcf | Yes | | |

## PRS Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | prs_test_code | Yes | Test code that defines config files and assay | |
| Array[Array[File]] | model_manifests | Yes | Adjustment model manifest files from MakemodelModelWorkflow WDL | |
| File | conditions_config | Yes | TSV file with condition-specific data on bins, odds-ratios, codes, etc. | |
| File | percentiles | Yes | TSV file with PRS scores corresponding to percentiles (0-100) for each condition | |
| File | score_weights | No | Weights for PRS scores in order to create PRSmix scores; Required to only if mix_before_adjustment is false |
| Boolean | mix_before_adjustment | No | If `true`, mix raw scores before adjusting them |  |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| File | renaming_lookup | No | Mapping file for renaming chromosomes | "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv" |
| String | ubuntu_docker_image | No | Ubuntu Docker image | "ubuntu:latest" |
| String | prs_docker_image | No | Docker image with PRS scripts | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | glimpse_vcf | Always | Output imputed VCF from running GLIMPSE |
| File | glimpse_qc_metrics | Always | QC metrics from running GLIMPSE |
| Array[File] | prs_score | Always | PRS mixed and adjusted scores |
| File | risk_summary | Always | TSV files with risk percentile, group, and odds ratio for the sample |
