version 1.0

import "VcfToBed.wdl"

workflow PCAWorkflow {
  input {
    File vcf
    File variants
    Int  memory   = 8
    Int  nthreads = 16
  }

  call VcfToBed.VcfToBed {
    input:
        vcf      = vcf
      , variants = variants
  }

  call BedToPca {
    input:
        bed      = VcfToBed.bed
      , bim      = VcfToBed.bim
      , fam      = VcfToBed.fam
  }

  output {
    File components   = BedToPca.components
    File eigenvalues  = BedToPca.eigenvalues
    File eigenvectors = BedToPca.eigenvectors
    File loadings     = BedToPca.loadings
    File meansd       = BedToPca.meansd
    File variance     = BedToPca.variance
  }
}


task BedToPca {
  input {
    File   bed
    File   bim
    File   fam
    Int    ndimensions = 20
    Int    nthreads    = 16

    Int    memory      = 8
    String docker      = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
  }

  Int    storage  = 12 * ceil(size(bed, "GB")) + 20

  String workdir  = "work"
  String basename = "pca"
  String prefix   = workdir + "/" + basename

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace

  mkdir --parents '~{workdir}'
  cp '~{bed}' '~{prefix}.bed'
  cp '~{bim}' '~{prefix}.bim'
  cp '~{fam}' '~{prefix}.fam'

  ~/flashpca/flashpca                      \
      --bfile      '~{prefix}'             \
      --numthreads ~{nthreads}             \
      --ndim       ~{ndimensions}          \
      --outpc     '~{prefix}'.components   \
      --outval    '~{prefix}'.eigenvalues  \
      --outvec    '~{prefix}'.eigenvectors \
      --outload   '~{prefix}'.loadings     \
      --outmeansd '~{prefix}'.meansd       \
      --outpve    '~{prefix}'.variance
  >>>

  output {
    File components   = prefix + ".components"
    File eigenvalues  = prefix + ".eigenvalues"
    File eigenvectors = prefix + ".eigenvectors"
    File loadings     = prefix + ".loadings"
    File meansd       = prefix + ".meansd"
    File variance     = prefix + ".variance"
  }

  runtime {
    docker: docker
    disks:  "local-disk " + storage + " HDD"
    memory: memory + " GB"
  }
}
