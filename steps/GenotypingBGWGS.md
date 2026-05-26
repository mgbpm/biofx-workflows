# GenotypingBGWGSWorkflow

This workflow performs variant genotyping for genome samples, producing a fully annotated VCF file. It optionally fetches the input CRAM from a GCP data location before running GATK HaplotypeCaller and GenotypeGVCFs.

## Input Parameters

| Type | Name | Req'd | Description |
| :--- | :--- | :---: | :---------- |
| String | sample_id | Yes | Full sample identifier |
| String | accession_id | Yes | Sample accession (PM#) |
| String | gcp_project_id | No | GCP project ID (default: "mgb-lmm-gcp-infrast-1651079146") |
| String | data_location | No | Data location for fetching CRAM files |
| Array[String] | fetch_cram_filter_keys | No | Keys used to filter CRAM files during fetch (default: [accession_id, sample_id]) |
| Boolean | fetch_cram | No | Whether to fetch the CRAM from a remote location (default: false) |
| Array[FileMatcher] | fetch_cram_file_matchers | No | File matchers for CRAM fetch |
| Boolean | fetch_files_verbose | No | Enable verbose logging for file fetching (default: false) |
| String | workspace_name | Yes | Workspace name for file fetching |
| Int | fetch_disk_size | No | Disk size in GB for fetch task (default: 75) |
| String | orchutils_docker_image | No | Docker image for orchestration utils (default: "us-docker.pkg.dev/mgbpmbiofx/orchestration-utils/orchestration-utils:20250203") |
| File | input_cram | No | Sample CRAM/BAM file (optional if fetch_cram is true) |
| File | input_crai | No | Corresponding CRAI/BAI file (optional if fetch_cram is true) |
| String | test_code | Yes | Workflow test-code |
| File | reference_fasta | Yes | HG38 Reference FASTA file |
| File | reference_fasta_fai | Yes | HG38 Reference FASTA index file |
| File | reference_dict | Yes | HG38 Reference dict file |
| File | roi_bed | Yes | Bed file with the workflow region of interest |
| File | dbsnp | Yes | HG38 DBSNP VCF file |
| File | dbsnp_vcf_index | Yes | HG38 DBSNP VCF index file |
| String | gatk_path | No | Path to GATK executable (default: /gatk/gatk) |
| String | genotyping_docker_image | Yes | Name of docker image for the workflow |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :---------- |
| File | annotated_vcf_file | Always | VCF file with genotypes and caller annotation |

## Workflow Steps

1. **FetchCram** *(conditional: only when fetch_cram is true)*  
   Fetches CRAM/CRAI files from a GCP data location using the `FetchFilesTask` from `FileUtils.wdl`.

2. **HaplotypeCallerTask**  
   Runs GATK HaplotypeCaller to generate a GVCF (BP_RESOLUTION mode), then runs GATK GenotypeGVCFs to produce an all-calls VCF.

3. **AddAnnotationsTask**  
   Annotates the VCF with caller/version information using `annotate_with_caller.py`.

## Example Output File Names

- `<accession_id>_<sample_id>_<test_code>.g.vcf`
- `<accession_id>_<sample_id>_<test_code>.allcalls.vcf`
- `<accession_id>_<sample_id>_<test_code>.annotated.vcf`

---