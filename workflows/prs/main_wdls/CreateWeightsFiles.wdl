version 1.0

workflow CreatePrsWeightsWorkflow {
    input {
        String pgs_id
        File b38_lookup
        File b37_lookup
        File chain_file
        String prs_docker_image
  }

    call DownloadWeightsFileTask {
        input:
            pgs_id = pgs_id,
            docker_image = prs_docker_image
    }

    if (DownloadWeightsFileTask.build == "37") {
        call CreateWeightsTask as CreateB37Weights {
            input:
                pgs_id = pgs_id,
                build = "37",
                input_weights = select_first([DownloadWeightsFileTask.b37_weights_file]),
                lookup_file = b37_lookup,
                docker_image = prs_docker_image
        }

        call LiftoverWeightsTask {
            input:
                pgs_id = pgs_id,
                input_weights = CreateB37Weights.output_weights,
                bed_file = select_first([CreateB37Weights.b37_bed_file]),
                chain_file = chain_file,
                docker_image = prs_docker_image
        }
    }

    if (DownloadWeightsFileTask.build == "38") {
        call CreateWeightsTask as CreateB38Weights {
            input:
                pgs_id = pgs_id,
                build = "38",
                input_weights = select_first([DownloadWeightsFileTask.b38_weights_file]),
                lookup_file = b38_lookup,
                docker_image = prs_docker_image
        }
    }

    output {
        File downloaded_weights = select_first([DownloadWeightsFileTask.b37_weights_file, DownloadWeightsFileTask.b38_weights_file])
        File output_weights = select_first([LiftoverWeightsTask.output_weights, CreateB38Weights.output_weights])
    }
}

task DownloadWeightsFileTask {
    input {
        String pgs_id
        String docker_image
        Int disk_size = 10
        Int mem_size = 2
        Int preemptible = 1
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
        File? b38_weights_file = "~{pgs_id}_hmPOS_GRCh38.txt.gz"
        File? b37_weights_file = "~{pgs_id}_hmPOS_GRCh37.txt.gz"
        String build = read_string("build.txt")
    }
}

task CreateWeightsTask {
    input {
        String pgs_id
        String build
        File input_weights
        File lookup_file
        String docker_image
        Int addldisk = 15
        Int bootdisk_size = 20
        Int mem_size = 4
        Int preemptible = 1
    }

    Int file_size = ceil(size(lookup_file, "GB") + size(input_weights, "GB"))
    Int final_disk_size = file_size + addldisk

    command <<<
        set -euxo pipefail

        mkdir -p OUTPUT

        /mgbpmbiofx/packages/biofx-prs/create_weights/make_weights_file.py \
                "~{input_weights}" \
                "~{lookup_file}" \
                "OUTPUT"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        bootDiskSizeGb: bootdisk_size
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_weights = "~{pgs_id}_hmPOS_GRCh~{build}.weights.tsv"
        File? b37_bed_file = "~{pgs_id}_hmPOS_GRCh~{build}.original.bed"
        File errors = "~{pgs_id}_hmPOS_GRCh~{build}.errors.tsv"
    }
}

task LiftoverWeightsTask {
    input {
        String pgs_id
        File input_weights
        File bed_file
        File chain_file
        String docker_image
        Int addldisk = 10
        Int mem_size = 4
        Int preemptible = 1
    }

    Int file_size = ceil(size(chain_file, "GB") + size(input_weights, "GB") + size(bed_file, "GB"))
    Int final_disk_size = file_size + addldisk

    command <<<
        # Liftover GRCh37 bed to GRCh38
        /mgbpmbiofx/packages/biofx-prs/create_weights/liftOver \
            "~{bed_file}" \
            "~{chain_file}" \
            "lifted.bed" \
            "unmapped.bed"

        # Liftover weights using GRch38 bed
        /mgbpmbiofx/packages/biofx-prs/create_weights/liftover_weights.py \
            "lifted.bed" \
            "unmapped.bed" \
            "~{input_weights}" \
            "."
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_weights = "~{pgs_id}_hmPOS_GRCh38.weights.tsv"
    }
}