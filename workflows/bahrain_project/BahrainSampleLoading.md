# Bahrain Sample Loading Workflow

The Bahrain Sample Loading Workflow is a continuation of running the monogenic or screening pipeline after the Bahrain Sample Prep Workflow. It will annotate and load individual sample VCFs to FAST.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | VCF file for annotation and loading to FAST |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20240625" |
| String | qceval_project_type | No | The type of rules to apply for the QC evaluation task, one of "WGS", "WGS_DRAGEN", "WES" or "NONE" | "WGS" |
| String | qceval_docker_image | No | The name of the Docker image to run the QC evaluation task | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20230511" |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| Boolean | has_haploid_sites | No | If true, modify the VCF file headers prior to FAST load to work around lack of support Number=G fields and haploid sites | false |
| String | sample_data_load_config_name | No | The FAST load configuration name for the sample data VCF | "Sample_VCF_PPM_Eval" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | qceval_vcf_gz | Always | Target VCF file annotated with QC Evaluation |
| String | fast_sample_data_name | Always | The name used when loading the sample to FAST |
| String | qceval_data_load_id | Always | The ID associated with loading the QCEval vcf to FAST |
