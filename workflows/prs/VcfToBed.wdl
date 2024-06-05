version 1.0

task VcfToBed {
  input {
    File   vcf
    File   variants
    File?  subset
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

  if '~{if defined(subset) then "true" else "false"}'
  then
      EXTRACT=( --extract-intersect '~{variants}' '~{subset}' )
  else
      EXTRACT=( --extract-intersect '~{variants}' )
  fi

  /plink2                           \
      --allow-extra-chr             \
      ${EXTRACT[@]+"${EXTRACT[@]}"} \
      --make-bed                    \
      --new-id-max-allele-len       \
          1000                      \
          missing                   \
      --out                         \
          '~{prefix}'               \
      --rm-dup                      \
          force-first               \
      --set-all-var-ids             \
          '@:#:$1:$2'               \
      --vcf                         \
          '~{vcf}'
  >>>

  output {
    File           bed = "~{prefix}.bed"
    File           bim = "~{prefix}.bim"
    File           fam = "~{prefix}.fam"
    Array[File] INPUTS = inputs
  }

  runtime {
    docker: docker
    disks:  "local-disk " + storage + " HDD"
    memory: memory + " GB"
  }
}
