version 1.0

# VariantListIgvScreenshots.wdl
#
# For each unique Biosample_ID in the input variant list, this workflow:
#
#   1. Looks up every matching source directory from the manifest file.
#   2. Copies the CRAM (and CRAI, MD5) files for that biosample from Wasabi S3
#      to a GCS staging bucket using CopySampleFilesWorkflow.
#   3. Generates an HTML IGV screenshot report—one per CRAM—covering all
#      variants listed for that biosample.
#
# All biosamples are processed in parallel (one scatter branch per biosample).
#
# Inputs
# ------
#   variant_list_tsv  : TSV with columns Variant_ID_VCF, Biosample_ID
#                       Variant_ID_VCF format: CHR-POS-REF-ALT  (e.g. "X-31126647-T-A")
#
#   manifest_tsv      : TSV with columns subject_id, path, filename, size
#                       subject_id is matched against Biosample_ID (exact, case-sensitive)
#                       path is a directory path relative to s3_prefix

import "../../steps/FileUtils.wdl" as FilesUtilsWf

workflow VariantListIgvScreenshots {

    input {
        # -----------------------------------------------------------------------
        # Input data files
        # -----------------------------------------------------------------------
        File   variant_list_tsv    # Variant_ID_VCF, Biosample_ID
        File   manifest_tsv        # subject_id, path, filename, size

        # -----------------------------------------------------------------------
        # Reference genome (required by igv-reports)
        # -----------------------------------------------------------------------
        File   ref_fasta
        File   ref_fasta_index     # Accompanying .fai index

        # -----------------------------------------------------------------------
        # IGV screenshot settings
        # -----------------------------------------------------------------------
        # Bases to display on either side of each variant region
        Int    igv_flanking = 50
        # Disk size (GB) for the IGV task — size for large CRAM localization
        Int    igv_disk_gb  = 200

        # -----------------------------------------------------------------------
        # Wasabi / S3 source configuration
        # source_location = s3_prefix + "/" + manifest.path
        # -----------------------------------------------------------------------
        String s3_prefix

        # -----------------------------------------------------------------------
        # CopySampleFilesWorkflow configuration
        # -----------------------------------------------------------------------
        String        staging_bucket
        String        gcp_project_id
        String        workspace_name
        # File extensions to copy; includes crai so the IGV task can find the index
        Array[String] cram_file_types   = ["cram", "crai"]

        # -----------------------------------------------------------------------
        # Docker images
        # -----------------------------------------------------------------------
        String igvreport_docker_image                 = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511"
        String orchutils_docker_image                 = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
        String variantlistigvscreenshots_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/wasabi-igv-screenshot:202608142"

        Int    preemptible = 1
    }

    # -------------------------------------------------------------------------
    # Step 1 — Parse inputs
    # Produces one variant TSV file and one source-paths file per unique Biosample_ID.
    # -------------------------------------------------------------------------
    call PrepSampleDataTask {
        input:
            variant_list_tsv = variant_list_tsv,
            manifest_tsv     = manifest_tsv,
            s3_prefix        = s3_prefix,
            docker_image     = variantlistigvscreenshots_docker_image,
            preemptible      = preemptible
    }

    # -------------------------------------------------------------------------
    # Step 2 — One scatter branch per unique Biosample_ID
    # -------------------------------------------------------------------------
    scatter (i in range(length(PrepSampleDataTask.biosample_ids))) {

        String biosample_id = PrepSampleDataTask.biosample_ids[i]
        File   variant_tsv  = PrepSampleDataTask.variant_tsv_files[i]

        # Source S3 directory paths for this biosample (one per manifest row)
        Array[String] source_paths = read_lines(PrepSampleDataTask.source_paths_files[i])

        # -------------------------------------------------------------------
        # Step 2a — Copy CRAM(s) from each matching manifest directory to GCS
        # One CopySampleFilesWorkflow call per source directory (inner scatter).
        #
        # Collision protection: each call receives a unique target subdirectory
        # keyed by biosample_id and the source-path index (j).  This prevents
        # CRAMs with identical filenames from different source directories from
        # overwriting each other in the staging bucket when flatten = true.
        #   e.g.  .../cram_stage/10000001/0/sample.cram
        #         .../cram_stage/10000001/1/sample.cram   ← distinct, not overwritten
        # -------------------------------------------------------------------
        scatter (j in range(length(source_paths))) {
            call FilesUtilsWf.CopyFilesTask as CopyCram {
                input:
                    file_match_keys        = [biosample_id],
                    source_location        = source_paths[j],
                    flatten                = true,
                    recursive              = true,
                    target_location        = staging_bucket + "/" + biosample_id + "/" + j,
                    docker_image           = orchutils_docker_image,
                    gcp_project_id         = gcp_project_id,
                    workspace_name         = workspace_name,
                    file_types             = cram_file_types,
                    verbose                = true
            }
        }

        # CopyFilesTask.local_files is only populated when files are copied to
        # local disk. Here target_location is a GCS bucket, so consume
        # target_files and coerce to File to force localization in the next task.
        Array[File] copied_cram_files = flatten(CopyCram.target_files)

        # -------------------------------------------------------------------
        # Step 2b — Generate IGV screenshots
        # All copied files across every directory are localized by WDL as File
        # inputs; the task finds CRAMs/CRAIs and runs create_report for each
        # CRAM against every variant in the biosample's TSV file.
        # -------------------------------------------------------------------
        call IgvReportFromVariantTsvTask {
            input:
                all_localized_files = copied_cram_files,
            variant_tsv         = variant_tsv,
                biosample_id        = biosample_id,
                ref_fasta           = ref_fasta,
                ref_fasta_index     = ref_fasta_index,
                igv_flanking        = igv_flanking,
                disk_gb             = igv_disk_gb,
                docker_image        = igvreport_docker_image,
                preemptible         = preemptible
        }
    }

    output {
        # Outer array: one entry per biosample
        # Inner array: one HTML report per CRAM found for that biosample
        Array[Array[File]] igv_reports = IgvReportFromVariantTsvTask.igv_report_htmls
    }
}

# =============================================================================
# Task: PrepSampleDataTask
#
# Reads the variant list TSV and manifest TSV and produces, for each unique
# Biosample_ID:
#
#   variants/{index}.tsv  — per-biosample TSV file (with header):
#                           CHR, START, END, REF, ALT, Biosample_ID
#
#   paths/{index}.txt     — One S3 source_location path per line, built from
#                           s3_prefix + manifest.path for every manifest row
#                           whose subject_id matches the Biosample_ID.
#
# Output arrays are guaranteed to be aligned: biosample_ids[i] corresponds to
# variant_tsv_files[i] and source_paths_files[i].
# =============================================================================

task PrepSampleDataTask {

    input {
        File   variant_list_tsv
        File   manifest_tsv
        String s3_prefix
        String docker_image
        Int    preemptible = 1
    }

    command <<<
        set -euxo pipefail

        $MGBPMBIOFXPATH/biofx-igv-screenshot/bin/prep_sample_data.py \
            '~{variant_list_tsv}' \
            '~{manifest_tsv}'     \
            '~{s3_prefix}'
    >>>

    runtime {
        docker:      "~{docker_image}"
        disks:       "local-disk 10 HDD"
        preemptible: preemptible
        memory:      "2 GB"
    }

    output {
        # Parallel arrays — biosample_ids[i] corresponds to
        # variant_tsv_files[i] and source_paths_files[i]
        Array[String] biosample_ids      = read_lines("biosample_ids.txt")
        # Collect generated files directly from task outputs to avoid any
        # dependency on intermediate path-manifest files.
        Array[File]   variant_tsv_files  = glob("variants/*.tsv")
        Array[File]   source_paths_files = glob("paths/*.txt")
    }
}

# =============================================================================
# Task: IgvReportFromVariantTsvTask
#
# Given the flat list of all files copied by CopySampleFilesWorkflow (CRAMs,
# CRAIs, MD5s), this task:
#
#   1. Identifies every .cram file and its matching .crai index.
#   2. For each CRAM, calls igv-reports' create_report using the per-biosample
#      variant TSV file (CHR, START, END, REF, ALT, Biosample_ID).
#   3. Outputs one HTML report per CRAM.
#
# CRAI matching logic
# -------------------
# Given a CRAM at any local path with basename "sample.cram", the task searches
# all_localized_files for a file whose path ends with either:
#   - /sample.cram.crai   (index co-named with CRAM, e.g. Dragen default)
#   - /sample.crai        (index with only the sample stem)
# =============================================================================

task IgvReportFromVariantTsvTask {

    input {
        # All files returned by CopySampleFilesWorkflow (CRAMs, CRAIs, MD5s)
        Array[File] all_localized_files
        # Per-biosample variant TSV produced by PrepSampleDataTask
        File        variant_tsv
        String      biosample_id
        File        ref_fasta
        File        ref_fasta_index
        Int         igv_flanking = 50
        Int         disk_gb      = 200
        String      docker_image
        Int         preemptible  = 1
    }

    command <<<
        set -euxo pipefail

        # Link reference FASTA index alongside the FASTA for igv-reports
        [ ! -f "~{ref_fasta}.fai" ] && ln -sf "~{ref_fasta_index}" "~{ref_fasta}.fai"

        # Write all localized file paths to a temp file for safe bash processing
        # (avoids word-splitting on paths that may contain spaces)
        ALL_FILES_LIST=~{write_lines(all_localized_files)}

        # Separate CRAMs and CRAIs from the mixed file list
        grep -i '\.cram$' "$ALL_FILES_LIST" | sort > cram_files.txt || true
        grep -i '\.crai$' "$ALL_FILES_LIST" | sort > crai_files.txt || true

        num_crams=$(wc -l < cram_files.txt)

        if [ "${num_crams}" -eq 0 ]; then
            echo "WARNING: No CRAM files found for biosample ~{biosample_id}." >&2
            touch "~{biosample_id}_no_crams.igvreport.html"
            exit 0
        fi

        # Count variant rows after removing header line
        num_variants=$(tail -n +2 "~{variant_tsv}" | wc -l)
        if [ "${num_variants}" -eq 0 ]; then
            echo "WARNING: Variant TSV file is empty for biosample ~{biosample_id}." >&2
            touch "~{biosample_id}_no_variants.igvreport.html"
            exit 0
        fi

        # Create a working directory for symlinks so igv-reports can find each
        # CRAM's index without needing write access to the localized file paths
        mkdir -p working

        while IFS= read -r cram_path; do
            cram_base=$(basename "${cram_path}")     # e.g. sample.cram
            cram_stem="${cram_base%.cram}"           # e.g. sample

            # Search for the CRAI that belongs to this CRAM.
            # Dragen convention: <name>.cram.crai; also accept <name>.crai
            crai_path=$(grep -E "(/${cram_base}\.crai$|/${cram_stem}\.crai$)" \
                        "$ALL_FILES_LIST" | head -1 || true)

            if [ -z "${crai_path}" ]; then
                echo "WARNING: No CRAI found for ${cram_path}; skipping." >&2
                continue
            fi

            # Symlink CRAM and its index into the working directory
            # igv-reports expects the index at <cram_path>.crai
            ln -sf "${cram_path}" "working/${cram_base}"
            ln -sf "${crai_path}" "working/${cram_base}.crai"

            out_html="~{biosample_id}_${cram_stem}.igvreport.html"

            create_report "~{variant_tsv}" "~{ref_fasta}"          \
                --sequence 1                                       \
                --begin    2                                       \
                --end      3                                       \
                --flanking ~{igv_flanking}                         \
                --info-columns CHR STAR END REF ALT Biosample_ID   \
                --tracks   "working/${cram_base}"                  \
                --output   "${out_html}"

        done < cram_files.txt
    >>>

    runtime {
        docker:      "~{docker_image}"
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
        memory:      "8 GB"
    }

    output {
        # One HTML report per CRAM; placeholder file if no CRAMs were found
        Array[File] igv_report_htmls = glob("*.igvreport.html")
    }
}
