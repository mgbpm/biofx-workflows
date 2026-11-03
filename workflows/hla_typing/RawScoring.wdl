version 1.0

workflow RawScoreWorkflow {
    input {
        File   input_vcf
        File   weights
        String plink2_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    }

    call ScoreVcf {
        input:
            vcf          = input_vcf,
            weights      = weights,
            docker_image = plink2_docker_image
    }

    output {
        File prs_score    = ScoreVcf.score
        File sites_scored = ScoreVcf.sites_scored
    }
}

task ScoreVcf {
    input {
        File    vcf
        File    weights
        String  basename = sub(basename(vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
        String  docker_image
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

        /plink2                                     \
            --vcf                   '~{vcf}'        \
            --memory                ~{plink_memory} \
            --new-id-max-allele-len 1000 missing    \
            --rm-dup exclude-all                    \
            --set-all-var-ids       '@:#:$r:$a'     \
            --write-snplist         allow-dups      \
            --out "~{basename}"            

        COLUMNS='maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums'
        /plink2                                  \
            --vcf "~{vcf}"                       \
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
            --out "~{basename}"
            
    >>>

    output {
        File score        = "~{basename}.sscore"
        File sites_scored = "~{basename}.sscore.vars"
    }

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{final_disk_size} HDD"
        memory: "~{runtime_memory} GB"
        preemptible: preemptible
    }
}