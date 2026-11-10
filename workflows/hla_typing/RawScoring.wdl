version 1.0

workflow RawScoreWorkflow {
    input {
        File?  input_vcf
        File?  input_pgen
        File?  input_pvar
        File?  input_psam
        File   weights
        String plink2_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    }

    call ScoreVcf {
        input:
            vcf     = input_vcf,
            pgen    = input_pgen,
            pvar    = input_pvar,
            psam    = input_psam,
            weights = weights,
            docker_image = plink2_docker_image
    }

    output {
        File prs_score        = ScoreVcf.score
        File sites_scored     = ScoreVcf.sites_scored
        File sites_not_scored = ScoreVcf.sites_not_scored
    }
}

task ScoreVcf {
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
            INPUT_FILES=(--vcf "${VCF}" dosage=DS)
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