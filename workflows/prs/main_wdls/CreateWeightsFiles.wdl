version 1.0

workflow CreatePrsWeightsWorkflow {
    input {
        Array[String] pgs_ids
        File b38_lookup
        File b37_lookup
        File chain_file
        String prs_docker_image
  }

    scatter (id in pgs_ids) {
        call CreateWeightsTask {
            input:
                pgs_id = id,
                b38_lookup = b38_lookup,
                b37_lookup = b37_lookup,
                chain_file = chain_file,
                docker_image = prs_docker_image
        }
    }

    output {
        Array[File] output_weights_files = CreateWeightsTask.weights_file
    }
}


task CreateWeightsTask {
    input {
        String pgs_id
        File b38_lookup
        File b37_lookup
        File chain_file
        String docker_image
        Int addldisk = 15
        Int mem_size = 2
        Int preemptible = 1
    }

    Int lookup_size = ceil(size(b37_lookup, "GB") + size(b38_lookup, "GB"))
    Int final_disk_size = lookup_size + addldisk

    command <<<
        set -euxo pipefail

        mkdir WORK
        mkdir OUTPUT

        # Download GRCh38 harmonized PGS Catalog weights file
        wget --secure-protocol=TLSv1_2 \
            -P WORK \
            "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/~{pgs_id}/ScoringFiles/Harmonized/~{pgs_id}_hmPOS_GRCh38.txt.gz"

        # If "liftover" is the source of lift over in the GRCh38 file, use the GRCh37 file instead
        if [[ $(zgrep -c "liftover" "WORK/~{pgs_id}_hmPOS_GRCh38.txt.gz") -ge 1 ]]; then
            rm "WORK/~{pgs_id}_hmPOS_GRCh38.txt.gz"

            wget --secure-protocol=TLSv1_2 \
                -P WORK \
                "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/~{pgs_id}/ScoringFiles/Harmonized/~{pgs_id}_hmPOS_GRCh37.txt.gz"
        fi

        # Make new weights file from PGS Catalog file
        if [ -f "WORK/~{pgs_id}_hmPOS_GRCh37.txt.gz" ]; then
            make_weights_file.py \
                "WORK/~{pgs_id}_hmPOS_GRCh37.txt.gz" \
                "~{b37_lookup}" \
                "WORK"

            # Liftover GRCh37 bed to GRCh38
            liftOver \
                "WORK/~{pgs_id}_hmPOS_GRCh37.original.bed" \
                "~{chain_file}" \
                "WORK/~{pgs_id}_hmPOS_GRCh38.lifted.bed" \
                "WORK/~{pgs_id}_hmPOS_GRCh37.unmapped.bed"

            # Liftover weights using GRch38 bed
            liftover_weights.py \
                "WORK/~{pgs_id}_hmPOS_GRCh38.lifted.bed" \
                "WORK/~{pgs_id}_hmPOS_GRCh37.unmapped.bed" \
                "WORK/~{pgs_id}_hmPOS_GRCh37.weights.tsv" \
                "OUTPUT"
        else
            make_weights_file.py "WORK/~{pgs_id}_hmPOS_GRCh38.txt.gz" "~{b38_lookup}" "OUTPUT"
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File weights_file = "OUTPUT/~{pgs_id}_hmPOS_GRCh38.weights.tsv"
    }
}