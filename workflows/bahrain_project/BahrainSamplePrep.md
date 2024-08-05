# Bahrain Sample Prep Workflow

The Bahrain Sample Prep Workflow can be used to run the Bahrain monogenic and screening pipelines. Once either the monogenic or screening pipeline is specified, the workflow will filter Bahrain Genome Project data and create a joint VCF of the samples.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[String] | sample_ids | Yes | List of sample IDs, e.g. SM-MPPD5; are included in the vcf and cram file names in the iCGD Clinical workspace | |
| Array[String] | collaborator_sample_ids  | Yes | List of sample IDs, e.g. D-981108334-BH-4022-S1-A | |
| Array[String] | data_bucket | Yes | List of vcf (and cram for screening pipeline) data locations for each sample from iCGD Clinical workspace | |
| String | batch_name | Yes | Prefix for all FAST sample data names, e.g. BGP-BH-4022 or BGP-Batch10 | |
| File | target_roi_bed | Yes | BED file containing gene regions; the VCFs in the dataset will be filtered to these regions and a new merged vcf will be created | |
| String | pipeline_to_run | Yes | Either "monogenic" to run the monogenic pipeline or "screening" to run the screening pipeline | |
| String | python_docker_image | No | The name of the Python Docker image | "python:3.10" |
| String | bcftools_docker_image | No | The name of the bcftools Docker image for VCF manipulation | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17" |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20240625" |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| File | ref_fasta | Yes | The genome reference FASTA file | |
| File | ref_fasta_index | Yes | The genome reference FAST index file | |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | sample_vcfs | Always | Per sample VCFs after filtering and changing the header |
| Array[File] | sample_vcf_idx | Always | Per sample VCF index files that correspond to sample_vcfs |
| File | merged_vcf | Always | A merged VCF of sample_vcfs |
| File | collective_vcf | When running the screening pipeline |  A VCF with a single "fake" sample but contains all the variants in the merged VCF |
