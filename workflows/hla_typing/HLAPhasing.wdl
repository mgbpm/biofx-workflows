version 1.0

workflow HlaPhasingWorkflow {
    input {
        File   input_alleles
        String prs_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:dev-hla"
  }

    call PhasingTask {
        input:
            input_alleles = input_alleles,
            docker_image  = prs_docker_image
    }

    output {
        File scored_haplotypes = PhasingTask.scored_haplotypes
    }
}


task PhasingTask {
    input {
        File   input_alleles
        String output_basename = basename(input_alleles, ".tsv")
        String docker_image
        Int    disk_size = ceil(size(input_alleles, "GB")) + 10
        Int    mem_size = 2
        Int    preemptible = 1
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

        # Check for best pairs
        [[ -s "${PREFIX}.phasing.out_bestpairs" ]] || {
            echo "Error: no BESTPAIRS1 block found in ${PREFIX}.phasing.out" >&2
            exit 1
        }

        # Score haplotypes for T1DM
        Rscript ${DOCKERDIR}/bin/score_t1dm_maps.R \
            ${PREFIX}.phasing.out_bestpairs        \
            "~{output_basename}"                   \
            "OUTPUT"                               \
            ${DOCKERDIR}/etc/haplotype_weights.txt

        echo "Finished phasing and scoring for ~{output_basename}" >&2

        # ----------------------------------------------------------------
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File scored_haplotypes = "OUTPUT/~{output_basename}.phasing.out_scored"
        File phasing_log       = "OUTPUT/~{output_basename}.phasing.out"
    }
}

