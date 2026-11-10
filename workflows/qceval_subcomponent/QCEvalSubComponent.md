# QCEval Subcomponent Workflow

Runs the QC evaluation subcomponent of the BGWGS pipeline using WDL 1.0.
The workflow fetches a sample VCF from a storage location, filters it to the
target region of interest (ROI) BED file, then applies QCEval to the filtered VCF.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | gcp_project_id | No | GCP project used to retrieve secrets for file fetching | `mgb-lmm-gcp-infrast-1651079146` |
| String | workspace_name | Yes | Current Terra workspace name used for secret retrieval | |
| String | orchutils_docker_image | No | Docker image used to fetch sample files | `us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest` |
| String | bcftools_docker_image | No | Docker image used to filter the VCF | `us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17` |
| String | subject_id | Yes | Subject identifier used for file matching and output naming | |
| String | sample_id | Yes | Sample identifier used for file matching and output naming | |
| String | sample_data_location | Yes | Storage location from which to fetch the sample VCF | |
| Array[String] | fetch_vcf_filter_keys | No | Strings used to filter VCF file paths when custom file matchers are not defined | `[subject_id, sample_id]` |
| Array[FileMatcher]? | fetch_vcf_file_matchers | No | Custom file matching rules; when defined, replace the default VCF type and filter-key selection | |
| Boolean | fetch_files_verbose | No | Whether to enable verbose file-fetching output | `false` |
| Int | fetch_disk_size | No | Declared fetch disk size in GB; currently not passed to the fetch task | `75` |
| String | reference_build | No | Reference build label; currently not used by workflow tasks | `GRCh38` |
| File | ref_dict | Yes | Reference sequence dictionary; currently not used by workflow tasks | |
| File | ref_fasta | Yes | Reference FASTA passed to QCEval only for `BGE_DRAGEN_TP_BINNING` | |
| File | ref_fasta_index | Yes | Reference FASTA index passed to QCEval only for `BGE_DRAGEN_TP_BINNING` | |
| File | target_roi_bed | No | BED file defining the target ROI used to filter the VCF | `gs://lmm-reference-data/roi/targetROI_hg38_2023_08_24_withCHR.bed` |
| String | qceval_project_type | No | QCEval rules to apply. Supported values: `BGE_DRAGEN_TP_BINNING`, `WGS`, `WGS_DRAGEN`, `WES`, or `NONE` | `BGE_DRAGEN_TP_BINNING` |
| String | qceval_docker_image | No | Docker image used by the QCEval task | `us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20260521` |
| File? | thresholds | No | Thresholds TSV used for `BGE_DRAGEN_TP_BINNING` QCEval | `gs://lmm-reference-data/annotation/pmeval/thresholds_20250912.tsv` |
| File? | difficult_to_map_regions | No | Difficult-to-map regions archive used for `BGE_DRAGEN_TP_BINNING` QCEval | `gs://lmm-reference-data/annotation/pmeval/difficult_to_map_regions_20250912.tgz` |

For `BGE_DRAGEN_TP_BINNING`, both `thresholds` and `difficult_to_map_regions`
must be defined. The workflow explicitly fails if either is missing.
`ref_fasta` and `ref_fasta_index` are required workflow inputs for all project
types, but are passed to the QCEval task only for `BGE_DRAGEN_TP_BINNING`.

The `FileMatcher` struct is defined in [FileUtils.wdl](../../steps/FileUtils.wdl)
with fields `String file_type`, `Array[String]? path_matches`, and
`Boolean? required`.

## Output Parameters

| Type | Name | Description |
| :--- | :--- | :--- |
| File | target_vcf_gz | Compressed VCF filtered to the target ROI BED file, named `<subject_id>_<sample_id>.target.vcf.gz` |
| File | qceval_vcf_gz | Compressed, target-filtered VCF annotated by QCEval, named `<subject_id>_<sample_id>.qceval.vcf.gz` |

## Workflow Steps

1. Fetch files from `sample_data_location` without recursively searching
   subdirectories, using custom file matchers when defined or the default VCF
   type and filter keys otherwise.
2. Resolve the fetched VCF with `select_first([FetchVcf.vcf])`. The workflow
   fails if no VCF is returned.
3. Filter the fetched VCF to `target_roi_bed`.
4. For `BGE_DRAGEN_TP_BINNING`, check that the thresholds and regions archive
   are defined and supply the reference FASTA and index to QCEval.
5. Annotate the filtered VCF with QCEval using `qceval_project_type`.
