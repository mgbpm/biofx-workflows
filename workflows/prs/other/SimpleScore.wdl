version 1.0

workflow SimpleScore {
  input {
    File  vcf
    File  weights
    File? variants
  }

  call ScoreVcf {
    input:
      vcf      = vcf,
      weights  = weights,
      variants = variants
  }

  output {
    File scores   = ScoreVcf.scores
    File variants = ScoreVcf.variants
  }
}


task ScoreVcf {
  input {
    File   vcf
    File   weights
    File?  variants
    Int    memory   = 8
    String docker   = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
  }

  String outputdir = "final"
  String prefix    = outputdir + "/data"

  Int plink_memory = ceil(memory * 0.75 * 1000)
  Int storage      = 3 * ceil(size(vcf, "GB")) + 20

  String columns = "maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace

  # --------------------------------------------------------------------------

  CHROMOSOMES="$( mktemp )"

  tail --lines=+2 '~{weights}'     \
    | cut --fields=1 --delimiter=: \
    | sort --unique                \
    > "${CHROMOSOMES}"

  if grep --quiet --perl-regex '^chr' "${CHROMOSOMES}"
  then
      HEAD='chr'
  else
      HEAD=''
  fi

  if grep --quiet --perl-regex "^${HEAD}M$" "${CHROMOSOMES}"
  then
      ENCODING="${HEAD}M"
  else
      ENCODING="${HEAD}MT"
  fi

  rm "${CHROMOSOMES}"

  # --------------------------------------------------------------------------

  mkdir --verbose --parents "$( dirname '~{prefix}' )"

  if '~{if defined(variants) then "true" else "false"}'
  then
      EXTRACT=( --extract '~{variants}' )
  else
      EXTRACT=()
  fi

  /plink2                           \
      --allow-extra-chr             \
      ${EXTRACT[@]+"${EXTRACT[@]}"} \
      --new-id-max-allele-len       \
          1000                      \
          missing                   \
      --memory                      \
          ~{plink_memory}           \
      --out                         \
          '~{prefix}'               \
      --output-chr                  \
          "${ENCODING}"             \
      --score                       \
          '~{weights}'              \
          header                    \
          no-mean-imputation        \
          ignore-dup-ids            \
          list-variants             \
          cols='~{columns}'         \
      --set-all-var-ids             \
          '@:#:$r:$a'               \
      --vcf                         \
          '~{vcf}'
  >>>

  output {
    File scores   = "~{prefix}.sscore"
    File variants = "~{prefix}.sscore.vars"
    File log      = "~{prefix}.log"
  }

  runtime {
    docker: docker
    disks:  "local-disk " + storage + " HDD"
    memory: (memory + 2) + " GB"
  }
}
