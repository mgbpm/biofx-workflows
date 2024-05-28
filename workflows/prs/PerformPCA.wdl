version 1.0

workflow PCAWorkflow {
  input {
    File vcf
    File variants
    Int  memory   = 8
    Int  nthreads = 16
  }

  call VcfToBed {
    input:
        vcf      = vcf
      , variants = variants
  }

  call PerformPCA {
    input:
        bed      = VcfToBed.bed
      , bim      = VcfToBed.bim
      , fam      = VcfToBed.fam
  }

  output {
  }
}

task VcfToBed {
  input {
    File   vcf
    File   variants
    Int    memory   = 8
    String docker   = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
  }

  String prefix  = "data"
  Int    storage =  3 * ceil(size(vcf, "GB")) + 20

  String devdir     = 'DEV'
  String inputsdir  = devdir + '/INPUTS'
  String outputsdir = devdir + '/OUTPUTS'

  Array[String] inputs = [inputsdir + '/vcf',
                          inputsdir + '/variants']

  command <<<
  ### DEV START ###
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace

  mkdir --parents '~{inputsdir}' '~{outputsdir}'
  cp '~{vcf}'      "~{inputsdir}/vcf"
  cp '~{variants}' "~{inputsdir}/variants"

  # ------------------------------------------------------------------------
  ### DEV END ###

  /plink2                     \
      --allow-extra-chr       \
      --extract-intersect     \
          '~{variants}'       \
      --make-bed              \
      --new-id-max-allele-len \
          1000                \
          missing             \
      --out                   \
          '~{prefix}'         \
      --rm-dup                \
          force-first         \
      --set-all-var-ids       \
          '@:#:$1:$2'         \
      --vcf                   \
          '~{vcf}'
  >>>

  output {
    File bed = "~{prefix}.bed"
    File bim = "~{prefix}.bim"
    File fam = "~{prefix}.fam"
    Array[File] INPUTS = inputs
  }

  runtime {
    docker: docker
    disks: "local-disk " + storage + " HDD"
    memory: memory + " GB"
  }
}


task PerformPCA {
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
  String basename = "data"
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
      --outpc      '~{prefix}'.pc          \
      --outpve     '~{prefix}'.pc.variance \
      --outload    '~{prefix}'.pc.loadings \
      --outmeansd  '~{prefix}'.pc.meansd
  >>>

  output {
    File pc           = prefix  + ".pc"
    File variance     = prefix  + ".variance"
    File loadings     = prefix  + ".loadings"
    File meansd       = prefix  + ".meansd"
    File eigenvectors = workdir + "/eigenvectors.txt"
    File eigenvalues  = workdir + "/eigenvalues.txt"
  }

  runtime {
    docker: docker
    disks: "local-disk " + storage + " HDD"
    memory: memory + " GB"
  }
}
