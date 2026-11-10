version 1.0

task RawPrsTask {
    input {
        File?   vcf
        File?   pgen
        File?   pvar
        File?   psam
        File    weights
        String  basename     = if defined(vcf) then
                                   sub(basename(select_first([vcf])), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
                               else if defined(pgen) then
                                   sub(basename(select_first([pgen])), "\\.(pgen|PGEN)$", "")
                               else
                                   "placeholder"
        String  docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        Int     addldisk     = 20
        Int     mem_size     = 8
        Int     preemptible  = 1
    }

    Int base_memory     = if mem_size < 8 then 8 else mem_size
    Int plink_memory    = base_memory * 1000
    Int runtime_memory  = base_memory + 2
    Float dataset_size  = if defined(vcf) then size(select_first([vcf]), "GB")
                          else
                              (if defined(pgen) then size(select_first([pgen]), "GB") else 0.0) +
                              (if defined(pvar) then size(select_first([pvar]), "GB") else 0.0) +
                              (if defined(psam) then size(select_first([psam]), "GB") else 0.0)
    Int file_size       = 3 * ceil(dataset_size)
    Int final_disk_size = file_size + addldisk

    command <<<
        set -euxo pipefail

        VCF='~{if defined(vcf) then vcf else ""}'
        PGEN='~{if defined(pgen) then pgen else ""}'
        PVAR='~{if defined(pvar) then pvar else ""}'
        PSAM='~{if defined(psam) then psam else ""}'

        # Check input file compatibility
        if [[ -n "${VCF}" ]]; then
            if [[ -n "${PGEN}" || -n "${PVAR}" || -n "${PSAM}" ]]; then
                echo "Provide either a VCF or a PGEN/PVAR/PSAM dataset, not both." >&2
                exit 1
            fi
            # Check for dosage sum in both compressed and uncompressed VCFs
            if [[ "${VCF}" == *.gz || "${VCF}" == *.GZ || "${VCF}" == *.bgz || "${VCF}" == *.BGZ ]]; then
                VCF_GREP=zgrep
            else
                VCF_GREP=grep
            fi
            if "${VCF_GREP}" -q '^##FORMAT=<ID=DS[>,]' "${VCF}"; then
                INPUT_FILES=(--vcf "${VCF}" dosage=DS)
            else
                INPUT_FILES=(--vcf "${VCF}")
            fi
        elif [[ -n "${PGEN}" && -n "${PVAR}" && -n "${PSAM}" ]]; then
            INPUT_FILES=(--pgen "${PGEN}" --pvar "${PVAR}" --psam "${PSAM}")
        else
            echo "Provide either a VCF or a complete PGEN/PVAR/PSAM dataset." >&2
            exit 1
        fi

        # Extract variant IDs from the input dataset
        /plink2                                     \
            "${INPUT_FILES[@]}"                     \
            --memory                ~{plink_memory} \
            --new-id-max-allele-len 1000 missing    \
            --rm-dup exclude-all                    \
            --set-all-var-ids       '@:#:$r:$a'     \
            --write-snplist         allow-dups      \
            --allow-extra-chr                       \
            --out "~{basename}"

        # Score the input dataset
        COLUMNS='maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums'
        /plink2                                  \
            "${INPUT_FILES[@]}"                   \
            --extract "~{basename}.snplist"      \
            --score '~{weights}'                 \
                header                           \
                ignore-dup-ids                   \
                list-variants                    \
                no-mean-imputation               \
                cols="${COLUMNS}"                \
            --memory ~{plink_memory}             \
            --new-id-max-allele-len 1000 missing \
            --output-chr "MT"                    \
            --set-all-var-ids '@:#:$r:$a'        \
            --allow-extra-chr                    \
            --out "~{basename}"
        
        # Get a list of variants that were not scored
        comm -23                                 \
            <(sort -u "~{basename}.snplist")     \
            <(sort -u "~{basename}.sscore.vars") \
            > variants_unscored.txt

        # Keep the existing chromosome/position order for variant lists
        LC_ALL=C sort -t: -k1,1V -k2,2n \
            "~{basename}.sscore.vars" > "~{basename}.scored_vars.txt"
        LC_ALL=C sort -t: -k1,1V -k2,2n \
            variants_unscored.txt > "~{basename}.unscored_vars.txt"
    >>>

    output {
        File score                = "~{basename}.sscore"
        File sites_scored         = "~{basename}.scored_vars.txt"
        File sites_not_scored     = "~{basename}.unscored_vars.txt"
    }

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{final_disk_size} HDD"
        memory: "~{runtime_memory} GB"
        preemptible: preemptible
    }
}

task MixHlaPrsScoresTask {
    input {
        File   hla_score_file
        File   prs_score_file
        String docker_image = "ubuntu:latest"
        Int    addldisk = 10
        Int    mem_size = 2
    }

    Int input_size = ceil(size(hla_score_file, "GB") + size(prs_score_file, "GB"))
    Int disk_size  = input_size + addldisk

    command <<<
        set -euxo pipefail

        # Check that the sample sample IDs exist each input
        tail -n+2 "~{hla_score_file}" | cut -f1 | sort > "hla_sample_ids"
        tail -n+2 "~{prs_score_file}" | cut -f1 | sort > "prs_sample_ids"

        if ! cmp -s "hla_sample_ids" "prs_sample_ids"; then
            echo "Sample IDs in input HLA score file and PRS score file do not match." >&2
            exit 1
        fi

        printf -- "sample_id\thla_score\tprs_score\tfinal_score\n" > "t1dm_scores.tsv"

        while read line; do
            sample_id =$(echo "${line}" | cut -f1)

            hla_score=$(grep "${sample_id}" "~{hla_score_file}" | cut -f14)
            weighted_hla_score=$(awk -v x=0.79003 -v y=${hla_score} 'BEGIN {print x*y}')

            prs_score=$(grep "${sample_id}" "~{prs_score_file}" | cut -f4)
            weighted_prs_score=$(awk -v x=8.86562 -v y=${prs_score} 'BEGIN {print x*y}')

            merged_score=$(awk -v x=${weighted_hla_score} -v y=${weighted_prs_score} 'BEGIN {print x+y}')
            final_score=$(awk -v x=20 -v y=${merged_score} 'BEGIN {print x+y}')

            printf -- "${sample_id}\t${hla_score}\t${prs_score}\t${final_score}\n" \
                > "t1dm_scores.tsv"
        done < hla_sample_ids
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_size} GB"
    }

    output {
        File t1dm_scores = "t1dm_scores.tsv"
    }
}