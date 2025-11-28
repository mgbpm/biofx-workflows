version 1.0

# Forked from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

task MakePCAPlot {
  input {
    File   population_pcs
    File   target_pcs
    String basename
    String docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    Int    disk_size    = 100
    Int    mem_size     = 2
    Int    preemptible  = 1
  }

  command <<<
    Rscript -<< "EOF"
      library(dplyr)
      library(readr)
      library(ggplot2)

      population_pcs <- read_tsv("~{population_pcs}")
      target_pcs <- read_tsv("~{target_pcs}")

      ggplot(population_pcs, aes(x=PC1, y=PC2, color="Population")) +
        geom_point(size=0.1, alpha=0.1) +
        geom_point(data = target_pcs, aes(x=PC1, y=PC2, color="Target")) +
        labs(x="PC1", y="PC2") +
        theme_bw()

      ggsave(filename = "~{basename}_PCA_plot.png", dpi=300, width = 6, height = 6)

    EOF
  >>>

  output {
    File pca_plot = "~{basename}_PCA_plot.png"
  }

  runtime {
    docker: "~{docker_image}"
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task PerformPCA {
  input {
    File   bed
    File   bim
    File   fam
    String basename
    String docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
    Int    disk_size    = 400
    Int    mem_size     = 8
    Int    nthreads     = 16
    Int    preemptible  = 1
  }

  Int memory = 2 + (if mem_size < 8 then 8 else mem_size)

  command <<<
    WORKDIR="$( mktemp --directory )"
    PREFIX="${WORKDIR}/data"

    ln --symbolic '~{bed}' "${PREFIX}.bed"
    ln --symbolic '~{bim}' "${PREFIX}.bim"
    ln --symbolic '~{fam}' "${PREFIX}.fam"

    ~/flashpca/flashpca                      \
      --bfile      "${PREFIX}"               \
      --outpc      '~{basename}.pc'          \
      --outpve     '~{basename}.pc.variance' \
      --outload    '~{basename}.pc.loadings' \
      --outmeansd  '~{basename}.pc.meansd'   \
      --numthreads ~{nthreads}               \
      --ndim       20
  >>>

  output {
    File pcs          = "~{basename}.pc"
    File pc_loadings  = "~{basename}.pc.loadings"
    File mean_sd      = "~{basename}.pc.meansd"
    File pc_variance  = "~{basename}.pc.variance"
    File eigenvectors = "eigenvectors.txt"
    File eigenvalues  = "eigenvalues.txt"
  }

  runtime {
    docker: "~{docker_image}"
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{memory} GB"
    preemptible: preemptible
  }
}

task ProjectArray {
  input {
    File    bim
    File    bed
    File    fam
    File    pc_loadings
    File    pc_meansd
    String  basename
    Boolean orderIndependentCheck = false
    String  docker_image          = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
    Int     disk_size             = 400
    Int     mem_size              = 8
    Int     nthreads              = 16
    Int     preemptible           = 1
  }

  Int memory         = 2 + (if mem_size < 8 then 8 else mem_size)
  String postprocess = if orderIndependentCheck then "sort" else "cat"

  command <<<
    cp '~{bim}' '~{basename}.bim'
    cp '~{bed}' '~{basename}.bed'
    cp '~{fam}' '~{basename}.fam'

    cp '~{pc_loadings}' loadings.txt
    cp '~{pc_meansd}'   meansd.txt

    awk '{print $1}' loadings.txt | tail -n +2 > pcloadings_ids.txt
    awk '{print $1}' meansd.txt   | tail -n +2 > meansd_ids.txt

    diff pcloadings_ids.txt meansd_ids.txt > loadings_meansd_diff.txt
    rm pcloadings_ids.txt

    if [[ -s loadings_meansd_diff.txt ]]; then
        echo "PC loadings file and PC means file do not contain the same IDs (or in the same order); check your input files and run again." >&2
        exit 1
    fi
    rm loadings_meansd_diff.txt
    
    '~{postprocess}' meansd_ids.txt > pc_ids.txt
    rm meansd_ids.txt
    awk '{print $2}' '~{basename}.bim' | '~{postprocess}' > bim_ids.txt

    diff bim_ids.txt pc_ids.txt > bim_pc_diff.txt
    rm bim_ids.txt
    rm pc_ids.txt

    if [[ -s bim_pc_diff.txt ]]; then
        echo "IDs in .bim file are not the same as the IDs in the PCA files; check that you have the right files and run again." >&2
        exit 1
    fi
    rm bim_pc_diff.txt

    ~/flashpca/flashpca            \
      --verbose                    \
      --project                    \
      --numthreads ~{nthreads}     \
      --bfile      '~{basename}'   \
      --inmeansd   meansd.txt      \
      --inload     loadings.txt    \
      --outproj    "~{basename}_projections.txt"
  >>>

  output {
    File projections = "~{basename}_projections.txt"
  }

  runtime {
    docker: "~{docker_image}"
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{memory} GB"
    preemptible: preemptible
  }
}

task ArrayVcfToPlinkDataset {
  input {
    File   vcf
    File   pruning_sites
    File?  subset_to_sites
    String basename
    String docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    Int    addldisk     = 20
    Int    mem_size     = 8
    Int    preemptible  = 1
  }

  Int base_memory      = if mem_size < 8 then 8 else mem_size
  Int plink_mem        = base_memory * 1000
  Int runtime_memory   = base_memory + 2
  Int file_size        = 3 * ceil(size(vcf, "GB"))
  Int final_disk_size = addldisk + file_size
  String devdir        = 'DEV'
  String inputsdir     = devdir + '/INPUTS'
  String outputsdir    = devdir + '/OUTPUTS'
  Array[String] inputs = if defined(subset_to_sites)
                         then [inputsdir + '/vcf',
                               inputsdir + '/pruning_sites',
                               inputsdir + '/subset_to_sites']
                         else [inputsdir + '/vcf',
                               inputsdir + '/pruning_sites']
  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset
    set -o xtrace

    mkdir --parents '~{inputsdir}' '~{outputsdir}'
    cp '~{vcf}'           "~{inputsdir}/vcf"
    cp '~{pruning_sites}' "~{inputsdir}/pruning_sites"

    if '~{if defined(subset_to_sites) then "true" else "false"}'
    then
        cp '~{subset_to_sites}' "~{inputsdir}/subset_to_sites"
    fi

    /plink2 \
      --vcf ~{vcf} \
      --extract-intersect ~{pruning_sites} ~{subset_to_sites} \
      --allow-extra-chr \
      --set-all-var-ids '@:#:$r:$a' \
      --new-id-max-allele-len 1000 missing \
      --out ~{basename} \
      --make-bed \
      --memory ~{plink_mem} \
      --rm-dup force-first
  >>>

  output {
    File bed           = "~{basename}.bed"
    File bim           = "~{basename}.bim"
    File fam           = "~{basename}.fam"
    Array[File] INPUTS = inputs
  }

  runtime {
    docker: "~{docker_image}"
    disks: "local-disk ~{final_disk_size} HDD"
    memory: "~{runtime_memory} GB"
    preemptible: preemptible
  }
}
