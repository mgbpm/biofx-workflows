version 1.0

workflow CreatePrsWeightsWorkflow {
    input {
        String pgs_id
        File   b38_lookup_file
        File   b37_lookup_file
        File   chain_file
        String prs_docker_image
  }

    call DownloadWeightsFileTask {
        input:
            pgs_id = pgs_id,
            docker_image = prs_docker_image
    }

    if (DownloadWeightsFileTask.build == "37") {
        File b37_lookup = b37_lookup_file
    }
    if (DownloadWeightsFileTask.build == "38") {
        File b38_lookup = b38_lookup_file
    }

    call CreateWeightsTask {
        input:
            pgs_id = pgs_id,
            input_weights = DownloadWeightsFileTask.weights_file,
            lookup_file = select_first([b37_lookup, b38_lookup]),
            docker_image = prs_docker_image
    }

    if (DownloadWeightsFileTask.build == "37") {
        call LiftoverWeightsTask {
            input:
                pgs_id = pgs_id,
                input_weights = CreateWeightsTask.output_weights,
                bed_file = select_first([CreateWeightsTask.bed_file]),
                chain_file = chain_file,
                docker_image = prs_docker_image
        }
    }

    output {
        File output_weights = select_first([
            LiftoverWeightsTask.output_weights, CreateWeightsTask.output_weights
        ])
    }
}

task DownloadWeightsFileTask {
    input {
        String pgs_id
        String docker_image
        Int    disk_size   = 10
        Int    mem_size    = 2
        Int    preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # Download GRCh38 harmonized PGS Catalog weights file
        wget --secure-protocol=TLSv1_2 \
            "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/~{pgs_id}/ScoringFiles/Harmonized/~{pgs_id}_hmPOS_GRCh38.txt.gz"

        # If "liftover" is the source of lift over in the GRCh38 file, use the GRCh37 file instead
        if [[ $(zgrep -c "liftover" "~{pgs_id}_hmPOS_GRCh38.txt.gz") -ge 1 ]]; then
            rm "~{pgs_id}_hmPOS_GRCh38.txt.gz"

            wget --secure-protocol=TLSv1_2 \
                "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/~{pgs_id}/ScoringFiles/Harmonized/~{pgs_id}_hmPOS_GRCh37.txt.gz"
        
            printf -- "37" >> "build.txt"
        else
            printf -- "38" >> "build.txt"
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File   weights_file = select_first(["~{pgs_id}_hmPOS_GRCh37.txt.gz", "~{pgs_id}_hmPOS_GRCh38.txt.gz"])
        String build        = read_string("build.txt")
    }
}

task CreateWeightsTask {
    input {
        String pgs_id
        File   input_weights
        File   lookup_file
        String docker_image
        Int    addldisk    = 50
        Int    mem_size    = 16
        Int    preemptible = 1
    }

    Int file_size       = ceil(size(lookup_file, "GB") + size(input_weights, "GB"))
    Int final_disk_size = file_size + addldisk

    command <<<
        $PACKAGESDIR/biofx-prs/create_weights/make_weights_file.py \
            --scoring "~{input_weights}"                           \
            --lookup "~{lookup_file}"                              \
            --output-basename "~{pgs_id}"                          \
            --outputdir .
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File  output_weights = "~{pgs_id}.weights.tsv"
        File? bed_file       = "~{pgs_id}.original.bed"
        File  errors         = "~{pgs_id}.errors.tsv"
    }
}

task LiftoverWeightsTask {
    input {
        String pgs_id
        File   input_weights
        File   bed_file
        File   chain_file
        String docker_image
        Int    addldisk = 10
        Int    mem_size = 4
        Int    preemptible = 1
    }

    Int file_size       = ceil(size(chain_file, "GB") + size(input_weights, "GB") + size(bed_file, "GB"))
    Int final_disk_size = file_size + addldisk

    command <<<
        $PACKAGESDIR/biofx-prs/create_weights/liftOver \
            "~{bed_file}"                              \
            "~{chain_file}"                            \
            "lifted.bed"                               \
            "unmapped.bed"

        $PACKAGESDIR/biofx-prs/create_weights/liftover_weights.py \
            --lifted-bed "lifted.bed"                             \
            --liftover-errors "unmapped.bed"                      \
            --b37-weights "~{input_weights}"                      \
            --outputdir .
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_weights = "~{pgs_id}.weights.tsv"
    }
}