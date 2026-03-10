version 1.0

workflow FindVariantRegions {
    input {
        Array[File] weights
        File        pca_variants
        String      output_basename
        String      ubuntu_docker_image = "ubuntu:latest"
    }

    call FindRegionsTask {
        input:
            weights = weights,
            pca_variants = pca_variants,
            output_basename = output_basename,
            docker_image = ubuntu_docker_image
    }

    output {
        File regions = FindRegionsTask.regions
    }
}


task FindRegionsTask {
    input {
        Array[File] weights
        File        pca_variants
        String      output_basename
        String      docker_image
        Int         addldisk = 10
        Int         mem_size = 2
        Int         preemptible = 1
    }

    Int weights_size    = ceil(size(weights, "GB"))
    Int pca_size        = ceil(size(pca_variants, "GB"))
    Int final_disk_size = weights_size + pca_size + addldisk

    command <<<
        set -euxo pipefail

        cut -f1,2 -d ":" "~{pca_variants}" | sed 's/:/\t/g' > pca_variants.txt

        for file in '~{sep="' '" weights}'; do
            tail -n +2 ${file}     \
                | cut -f 1         \
                | cut -d ":" -f1,2 \
                | sed 's/^/chr/g'  \
                | sed 's/:/\t/g'   \
                >> weights_variants.txt
        done

        cat pca_variants.txt weights_variants.txt \
            | sort -u            \
            | sort -k1,1V -k2,2n \
            > "~{output_basename}.txt"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File regions = "~{output_basename}.txt"
    }
}