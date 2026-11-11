version 1.0

task FilterBam {
    input {
        File    input_file # can be CRAM or BAM
        File    input_file_idx
        File    ref_fasta       = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File    ref_fai         = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        File    ref_dict        = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
        File    intervals
        String? google_project
        String  output_basename = sub(basename(input_file), "\\.(cram|CRAM|bam|BAM)$", "")
        Int     cpu             = 2
        Int     num_threads     = 4
        Int     mem_size        = 4
        Int     addldisk        = 100
        Int     boot_disk_gb    = 10
        Int     preemptible     = 1
        String  docker_image    = "us.gcr.io/broad-gatk/gatk"
    }

    Int ref_size    = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB") + size(ref_dict, "GB"))
    Int sample_size = ceil(size(input_file, "GB") + size(input_file_idx, "GB"))
    Int disk_size   = addldisk + ref_size + ceil(size(input_file, "GB"))

    parameter_meta{
        intervals: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fai: {localization_optional: true}
        ref_dict: {localization_optional: true}
        input_file: {localization_optional: true}
        input_file_idx: {localization_optional: true}
    }

    command <<<
        set -euxo pipefail

        gatk PrintReads -R ~{ref_fasta} -I ~{input_file} -L ~{intervals} -O unsorted.bam \
            ~{if select_first([google_project, ""]) != "" then "--gcs-project-for-requester-pays " + select_first([google_project, ""]) else ""}

        gatk ValidateSamFile -I unsorted.bam

        samtools sort -@ ~{num_threads} unsorted.bam > "~{output_basename}.bam"

        gatk BuildBamIndex -I "~{output_basename}.bam"
    >>>

    runtime {
        docker:         docker_image
        bootDiskSizeGb: boot_disk_gb
        memory:         mem_size + " GB"
        disks:          "local-disk " + disk_size + " SSD"
        preemptible:    preemptible
        cpu:            cpu
    }

    output {
        File output_bam = "~{output_basename}.bam"
        File output_bai = "~{output_basename}.bai"
    }
}

task RunT1kTask {
    input {
        String input_format
        File   input_file
        File   input_index
        File   ref_fasta       = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File   ref_fai         = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        String preset_type     = "hla"
        String output_basename = sub(basename(input_file), "\\.(cram|CRAM|bam|BAM)$", "")
        String docker_image    = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/t1k:0.0.1"
        Int    mem_size        = ceil(16 + size(input_file, "GB") * 0.5)
        Int    addldisk        = 8
        Int    preemptible     = 1
    }

    Int sample_size        = ceil(size(input_file, "GB") + size(input_index, "GB"))
    Int ref_size           = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB"))
    Int storage_multiplier = if input_format == "cram" then 4 else 2
    Int final_disk_size    = (sample_size * storage_multiplier) + ref_size + addldisk

    command <<<
        set -euxo

        mkdir --parents "OUTPUT"
        mkdir --parents "TEMP"

        if [[ '~{input_format}' == 'cram' ]]
        then
            NTHREADS="$( nproc )"

            samtools view --bam              \
                --threads   "${NTHREADS}"    \
                --reference "~{ref_fasta}"   \
                --output    "TEMP/input.bam" \
                "~{input_file}"

            samtools index --bai           \
                --threads "${NTHREADS}"    \
                --output  "TEMP/input.bai" \
                "TEMP/input.bam"

            t1k "~{preset_type}" "TEMP/input.bam" "OUTPUT/~{output_basename}"
        else
            t1k "~{preset_type}" "~{input_file}" "OUTPUT/~{output_basename}"
        fi
    >>>

    runtime {
        docker:      docker_image
        disks:       "local-disk ~{final_disk_size} SSD"
        memory:      "~{mem_size}GB"
        preemptible: preemptible
    }

    output {
        File allele_tsv   = "OUTPUT/~{output_basename}_allele.tsv"
        File genotype_tsv = "OUTPUT/~{output_basename}_genotype.tsv"
        File allele_vcf   = "OUTPUT/~{output_basename}_allele.vcf"
    }
}

task PhasingTask {
    input {
        File   input_alleles
        String output_basename = basename(input_alleles, "_allele.tsv")
        String docker_image    = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:dev-hla"
        Int    disk_size       = ceil(size(input_alleles, "GB")) + 10
        Int    mem_size        = 2
        Int    preemptible     = 1
    }

    command <<<
        set -euxo pipefail

        mkdir -p OUTPUT

        PREFIX="OUTPUT/~{output_basename}"
        SAMPLE_ALLELE_CODES="${PREFIX}.allele_codes.txt"

        DOCKERDIR=${PACKAGESDIR}/biofx-prs/hla_imputation

        sites=(DRB1 DQA1 DQB1)

        # ----------------------------------------------------------------
        # Step 1: Extract allele codes
        # ----------------------------------------------------------------

        join_tab() { local IFS=$'\t'; printf '%s\n' "$*"; }

        # Look up an allele's code in the table; exit status 1 if not found.
        lookup_code() {
            awk -F'\t' -v a="$1" '
                { gsub(/\r/, ""); gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2) }
                $1 == a { print $2; found = 1; exit }
                END { exit !found }
            ' "${DOCKERDIR}/etc/allele_codes.tsv"
        }

        header=(ID)
        row=("~{output_basename}")

        for site in "${sites[@]}"; do
            alleles=()
            while IFS= read -r line; do
                alleles+=("$line")
            done < <(awk -v p="HLA-${site}*" 'index($1, p) == 1 { print $1 }' "~{input_alleles}" \
                     | cut -d: -f1,2 | tr '-' '_')

            if (( ${#alleles[@]} > 2 )); then
                echo "Warning: ${#alleles[@]} alleles for $site in ~{input_alleles}; using first two" >&2
            fi

            # Translate alleles to codes (unknown alleles get -1)
            codes=()
            if (( ${#alleles[@]} > 0 )); then
                for allele in "${alleles[@]:0:2}"; do
                    if code=$(lookup_code "$allele"); then
                        codes+=("$code")
                    else
                        echo "Warning: $allele not found in allele code table; writing -1 as allele" >&2
                        codes+=("-1")
                    fi
                done
            fi

            # Pad missing / single-allele sites with the dummy code
            while (( ${#codes[@]} < 2 )); do
                codes+=("-1")
            done

            header+=("${site}_1" "${site}_2")
            row+=("${codes[0]}" "${codes[1]}")
        done

        {
            join_tab "${header[@]}"
            join_tab "${row[@]}"
        } > "$SAMPLE_ALLELE_CODES"

        echo "Finished mapping allele codes" >&2

        # ----------------------------------------------------------------
        # Step 2: Phase alleles and score for T1DM
        # ----------------------------------------------------------------

        # Create phase input
        Rscript ${DOCKERDIR}/bin/create_phase_input.R \
            ${SAMPLE_ALLELE_CODES}                    \
            "~{output_basename}"                      \
            "OUTPUT"                                  \
            ${DOCKERDIR}/etc/known_haplotypes.txt

        # Check if PHASE input is missing/empty
        for f in "${PREFIX}.phasing.inp" "${PREFIX}.phasing.known"; do
            [[ -s $f ]] || { echo "Error: expected PHASE input $f is missing or empty" >&2; exit 1; }
        done

        # Run phasing
        ${DOCKERDIR}/src/PHASE         \
            -k${PREFIX}.phasing.known \
            -f1                       \
            -d1                       \
            -R0                       \
            ${PREFIX}.phasing.inp     \
            ${PREFIX}.phasing.out

        # Extract best haplotype pairs from PHASE output
        awk '/END BESTPAIRS1/ {f=0} f {print} /BEGIN BESTPAIRS1/ {f=1}' \
            "${PREFIX}.phasing.out"                                     \
            > "${PREFIX}.phasing.out_bestpairs"

        echo "Finished phasing for ~{output_basename}" >&2

        # ----------------------------------------------------------------
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File haplotype_pairs = "OUTPUT/~{output_basename}.phasing.out_bestpairs"
        File phasing_log       = "OUTPUT/~{output_basename}.phasing.out"
    }
}

task ScoreHaplotypesTask {
    input {
        File   haplotype_pairs
        String output_basename = basename(haplotype_pairs, ".phasing.out_bestpairs")
        String docker_image    = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:dev-hla"
        Int    disk_size       = ceil(size(haplotype_pairs, "GB")) + 10
        Int    mem_size        = 2
        Int    preemptible     = 1
    }

    command <<<
        set -euxo

        mkdir -p OUTPUT

        DOCKERDIR=${PACKAGESDIR}/biofx-prs/hla_imputation

        # Score haplotypes for T1DM
        Rscript ${DOCKERDIR}/bin/score_t1dm_maps.R \
            "~{haplotype_pairs}"              \
            "~{output_basename}"                   \
            "OUTPUT"                               \
            ${DOCKERDIR}/etc/haplotype_weights.txt

        echo "Finished scoring for ~{output_basename}" >&2
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File scored_haplotypes = "OUTPUT/~{output_basename}.phasing.out_scored"
    }
}

task CombineHaplotypeScoresTask {
    input {
        Array[File] scored_haplotype_files
        String      docker_image = "ubuntu:latest"
        Int         addldisk = 10
        Int         mem_size = 2
    }

    Int disk_size = ceil(size(scored_haplotype_files)) + addldisk

    command <<<
        set -euxo

        cat ~{sep=" " scored_haplotype_files} | sort -u > haplotype_scores.txt
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_size} GB"
    }

    output {
        File all_haplotype_scores = "haplotype_scores.txt"
    }
}