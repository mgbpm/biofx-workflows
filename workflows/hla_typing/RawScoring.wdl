version 1.0

workflow RawScoreWorkflow {
    input {
        File   input_vcf
        File   weights
    }

    call ScoreVcf {
        input:
            vcf          = input_vcf,
            weights      = weights
    }

    output {
        File prs_score        = ScoreVcf.score
        File sites_scored     = ScoreVcf.sites_scored
        File sites_not_scored = ScoreVcf.sites_not_scored
    }
}

task ScoreVcf {
    input {
        File    vcf
        File    weights
        String  basename     = sub(basename(vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
        String  docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        Int     addldisk     = 20
        Int     mem_size     = 8
        Int     preemptible  = 1
    }

    Int base_memory     = if mem_size < 8 then 8 else mem_size
    Int plink_memory    = base_memory * 1000
    Int runtime_memory  = base_memory + 2
    Int file_size       = 3 * ceil(size(vcf, "GB"))
    Int final_disk_size = file_size + addldisk

    command <<<
        set -euxo pipefail

        # Extract variant IDs from input VCF
        /plink2                                     \
            --vcf                   '~{vcf}'        \
            --memory                ~{plink_memory} \
            --new-id-max-allele-len 1000 missing    \
            --rm-dup exclude-all                    \
            --set-all-var-ids       '@:#:$r:$a'     \
            --write-snplist         allow-dups      \
            --allow-extra-chr                       \
            --out "~{basename}"

        # Score VCF
        COLUMNS='maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums'
        /plink2                                  \
            --vcf "~{vcf}" dosage=DS             \
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
        
        # Sort list of scored variants
        LC_ALL=C sort                 \
            -t: -k1,1V -k2,2n         \
            "~{basename}.sscore.vars" \
            > "~{basename}.scored_vars.txt"
        
        # Sort list of extracted variants from VCF
        LC_ALL=C sort             \
            -t: -k1,1V -k2,2n     \
            "~{basename}.snplist" \
            > "~{basename}.all_snps.txt"

        # Get list of variants that weren't scored
        comm -23                                 \
            <(sort -u "~{basename}.snplist")     \
            <(sort -u "~{basename}.sscore.vars") \
            > "~{basename}.not_scored"
        
        LC_ALL=C sort                \
            -t: -k1,1V -k2,2n        \
            "~{basename}.not_scored" \
            > "~{basename}.not_scored_vars.txt"

        # Get counts of variants
        wc -l < "~{basename}.not_scored_vars.txt" > not_scored.count
        wc -l < "~{basename}.scored_vars.txt" > scored.count
    >>>

    output {
        File score                = "~{basename}.sscore"
        File sites_scored         = "~{basename}.scored_vars.txt"
        Int  num_sites_scored     = read_int("scored.count")
        File sites_not_scored     = "~{basename}.not_scored_vars.txt"
        Int  num_sites_not_scored = read_int("not_scored.count")
    }

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{final_disk_size} HDD"
        memory: "~{runtime_memory} GB"
        preemptible: preemptible
    }
}