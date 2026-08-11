#!/usr/bin/env python3
"""
prep_sample_data.py

Parses a variant list TSV and a sample manifest TSV, then writes per-biosample
output files consumed by the VariantListIgvScreenshots WDL workflow.

Usage
-----
    prep_sample_data.py VARIANT_TSV MANIFEST_TSV S3_PREFIX

Arguments
---------
    VARIANT_TSV
        Tab-separated file with columns:
            Variant_ID_VCF   — CHR-POS-REF-ALT  (e.g. "X-31126647-T-A")
            Biosample_ID     — internal sample identifier (e.g. "10000001")
            Predicted_Impact — free-text impact annotation

    MANIFEST_TSV
        Tab-separated file with columns (at minimum):
            subject_id  — matched exactly against Biosample_ID
            path        — directory path relative to S3_PREFIX
            filename    — unused
            size        — unused

    S3_PREFIX
        S3 URI prefix prepended to each manifest path to build the
        source_location passed to CopySampleFilesWorkflow
        (e.g. "s3://prod-biobank-cram-2023-1").

Outputs  (written to the current working directory)
-------
    biosample_ids.txt        — sorted list of unique Biosample_IDs, one per line

    variants/{index}.bed     — 5-column BED file per biosample (with header):
                               chr, start (0-based), end, Predicted_Impact, Variant_ID_VCF

    paths/{index}.txt        — S3 source_location paths for that biosample,
                               one per line (empty if no manifest match found)

    variant_bed_paths.txt    — absolute paths to each variants/{index}.bed file,
                               in the same order as biosample_ids.txt
                               (consumed by WDL Array[File] output declaration)

    source_paths_paths.txt   — absolute paths to each paths/{index}.txt file,
                               in the same order as biosample_ids.txt
                               (consumed by WDL Array[File] output declaration)

Notes
-----
    - Chromosome labels are normalised to UCSC style: "X" -> "chrX", "MT"/"M" -> "chrM".
    - VCF positions (1-based) are converted to 0-based half-open BED intervals.
    - Output files use zero-padded numeric prefixes so that filesystem sort order
      matches the biosample_ids.txt ordering, which is required for WDL glob/
      read_lines alignment.
    - Input files may have a UTF-8 BOM; both are handled transparently.
"""

import csv
import os
import sys


def parse_variants(variant_tsv):
    """Return dict: biosample_id -> list of (chrom, bed_start, bed_end, impact, vid)."""
    variants_by_biosample = {}

    with open(variant_tsv, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            bid = row["Biosample_ID"].strip()
            vid = row["Variant_ID_VCF"].strip()
            impact = row["Predicted_Impact"].strip().replace("\t", " ")

            # Split on first 3 dashes to handle alleles that may contain '-'
            parts = vid.split("-", 3)
            if len(parts) < 4:
                print(f"WARNING: Cannot parse Variant_ID_VCF: {vid!r}", file=sys.stderr)
                continue

            chrom_raw, pos_str, ref, alt = parts

            # Normalise to UCSC chromosome names
            if chrom_raw.upper() in ("MT", "M"):
                chrom = "chrM"
            else:
                chrom = "chr" + chrom_raw

            try:
                pos = int(pos_str)
            except ValueError:
                print(f"WARNING: Invalid position in {vid!r}", file=sys.stderr)
                continue

            # Convert 1-based VCF position to 0-based half-open BED interval.
            # For SNVs: [pos-1, pos); IGV flanking handles the display window.
            bed_start = pos - 1
            bed_end = pos

            variants_by_biosample.setdefault(bid, []).append(
                (chrom, bed_start, bed_end, impact, vid)
            )

    return variants_by_biosample


def parse_manifest(manifest_tsv, known_biosample_ids, s3_prefix):
    """Return dict: biosample_id -> set of full S3 source_location paths."""
    paths_by_biosample = {}

    with open(manifest_tsv, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = row["subject_id"].strip()
            if sid not in known_biosample_ids:
                continue
            path = row["path"].strip()
            # Exactly one slash between prefix and path regardless of trailing/leading slashes
            full_path = s3_prefix + "/" + path.strip("/")
            paths_by_biosample.setdefault(sid, set()).add(full_path)

    return paths_by_biosample


def write_outputs(biosample_ids, variants_by_biosample, paths_by_biosample):
    """Write all per-biosample output files and the WDL path-manifest files."""
    os.makedirs("variants", exist_ok=True)
    os.makedirs("paths", exist_ok=True)

    bed_abs_paths = []
    paths_abs_paths = []

    for idx, bid in enumerate(biosample_ids):
        prefix = f"{idx:06d}"

        # BED file — header required by igv-reports for --info-columns
        bed_path = os.path.abspath(f"variants/{prefix}.bed")
        with open(bed_path, "w") as fh:
            fh.write("chr\tstart\tend\tPredicted_Impact\tVariant_ID_VCF\n")
            for chrom, s, e, impact, vid in variants_by_biosample[bid]:
                fh.write(f"{chrom}\t{s}\t{e}\t{impact}\t{vid}\n")
        bed_abs_paths.append(bed_path)

        # Source paths file — one S3 directory path per line (may be empty)
        paths_path = os.path.abspath(f"paths/{prefix}.txt")
        with open(paths_path, "w") as fh:
            for p in sorted(paths_by_biosample.get(bid, [])):
                fh.write(p + "\n")
        paths_abs_paths.append(paths_path)

    # WDL Array[File] output manifests — one absolute path per line,
    # aligned with biosample_ids.txt ordering
    with open("biosample_ids.txt", "w") as fh:
        fh.write("\n".join(biosample_ids) + "\n")

    with open("variant_bed_paths.txt", "w") as fh:
        fh.write("\n".join(bed_abs_paths) + "\n")

    with open("source_paths_paths.txt", "w") as fh:
        fh.write("\n".join(paths_abs_paths) + "\n")


def main():
    if len(sys.argv) != 4:
        print(
            f"Usage: {sys.argv[0]} VARIANT_TSV MANIFEST_TSV S3_PREFIX", file=sys.stderr
        )
        sys.exit(1)

    variant_tsv = sys.argv[1]
    manifest_tsv = sys.argv[2]
    s3_prefix = sys.argv[3].rstrip("/")

    variants_by_biosample = parse_variants(variant_tsv)

    if not variants_by_biosample:
        print(
            "ERROR: No variants could be parsed from variant list TSV.", file=sys.stderr
        )
        sys.exit(1)

    paths_by_biosample = parse_manifest(
        manifest_tsv, set(variants_by_biosample), s3_prefix
    )

    biosample_ids = sorted(variants_by_biosample.keys())

    write_outputs(biosample_ids, variants_by_biosample, paths_by_biosample)

    print(f"Prepared data for {len(biosample_ids)} biosample(s).", file=sys.stderr)

    unmatched = set(variants_by_biosample) - set(paths_by_biosample)
    if unmatched:
        print(
            f"WARNING: {len(unmatched)} biosample(s) had no matching manifest entry "
            f"and will produce no IGV screenshots: {sorted(unmatched)}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
