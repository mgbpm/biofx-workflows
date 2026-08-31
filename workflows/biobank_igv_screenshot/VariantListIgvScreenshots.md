## Variant List IGV Screenshots Workflow Module

`VariantListIgvScreenshots.wdl` generates IGV HTML reports from a variant list TSV and a sample manifest TSV.

For each unique `Biosample_ID` in the variant list, the workflow:

1. Builds one per-biosample variant TSV from `Variant_ID_VCF` values.
2. Finds all matching source directories from the manifest (`subject_id == Biosample_ID`).
3. Copies CRAM/CRAI files from S3/Wasabi paths into GCS staging locations.
4. Creates one IGV report per CRAM using `igv-reports`.

The workflow supports multiple manifest matches per biosample and emits distinct report filenames to avoid overwrite collisions when CRAM basenames are repeated.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | variant_list_tsv | Yes | Variant list TSV with columns `Variant_ID_VCF` and `Biosample_ID`. `Variant_ID_VCF` format: `CHR-POS-REF-ALT` | |
| File | manifest_tsv | Yes | Manifest TSV with columns `subject_id`, `path`, `filename`, `size` | |
| File | ref_fasta | Yes | Reference FASTA used by igv-reports | |
| File | ref_fasta_index | Yes | `.fai` index for `ref_fasta` | |
| Int | igv_flanking | No | Bases displayed on each side of variant interval in IGV report | 50 |
| Int | igv_disk_gb | No | Local disk size (GB) for report generation task | 200 |
| String | s3_prefix | Yes | Prefix prepended to manifest `path` to create source directories | |
| String | staging_bucket | Yes | GCS target prefix for copied CRAM/CRAI files | |
| String | gcp_project_id | Yes | Project ID used by file copy task | |
| String | workspace_name | Yes | Workspace name used by file copy task | |
| Array[String] | cram_file_types | No | File extensions copied from source directories | ["cram", "crai"] |
| String | igvreport_docker_image | No | Docker image for IGV report generation task | `us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511` |
| String | orchutils_docker_image | No | Docker image for copy task utilities | `us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest` |
| String | variantlistigvscreenshots_docker_image | No | Docker image containing `prep_sample_data.py` | `us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/wasabi-igv-screenshot:202608212` |
| Int | preemptible | No | Number of preemptible retries for tasks | 1 |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[Array[File]] | igv_reports | Always | Nested report outputs. Outer index is biosample, inner index is per-CRAM report file |
| Array[File] | all_igv_reports | Always | Flattened list of all report HTML files across all biosamples |

## Variant List Format

`variant_list_tsv` must have the following header columns:

- `Variant_ID_VCF`
- `Biosample_ID`

Example rows:

```text
Variant_ID_VCF	Biosample_ID
X-31126647-T-A	12345678
X-31126647-T-A	12345679
```

## Notes

- The prep task creates `variants/{index}.tsv` and `paths/{index}.txt` files aligned by biosample index.
- Report generation pairs each CRAM with a CRAI in the same directory when possible, then falls back to basename matching.
- Report filenames include a per-task index so repeated CRAM basenames generate multiple outputs instead of overwriting.
