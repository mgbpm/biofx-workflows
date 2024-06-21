version 1.0

import "PCATasks.wdl"
import "ScoringTasks.wdl"
import "TrainAncestryAdjustmentModel.wdl"
import "Structs.wdl"

struct RunSpec {
  File        vcf
  String      nvariants
  String      nsamples
  String      size
  Array[File] subset
  Int         memory
}

workflow Temp {
  input {
    Array[RunSpec] specs
  }

  scatter (spec in specs) {
    call VcfToBed {
      input:
        spec = spec
    }
  }
}

task VcfToBed {

  input {
    RunSpec spec
  }

  Int         storage   = 6 * ceil(size(vcf, "GB")) + 20
  File        vcf       = spec.vcf
  Int         memory    = spec.memory
  String      nvariants = spec.nvariants
  String      nsamples  = spec.nsamples
  String      size_     = spec.size
  String      results   = "results"
  Array[File] subset    = spec.subset
  Int         nsubset   = length(subset)

  command <<<
  apt-get update
  apt-get install --yes apt-utils
  apt-get install --yes time

  RESULTS='~{results}'

  # set -o errexit
  # set -o pipefail
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  (
      exec 1>>"${RESULTS}"
      printf -- '\n'
      printf -- 'NVARIANTS: %d; ' ~{nvariants}
      printf -- 'NSAMPLES: %d; '  ~{nsamples}
      printf -- 'SIZE: %d'        ~{size_}
      printf -- '\n\n'

      printf -- 'STORAGE: %d\n'   ~{storage}
      printf -- 'df BEFORE:\n'
      df --block-size=1 /cromwell_root
      df --human        /cromwell_root
  )

  if (( ~{nsubset} > 0 ))
  then
      EXTRACT=( --extract-intersect '~{sep="' '" subset}' )
  else
      EXTRACT=()
  fi

  PREFIX=data

  /usr/bin/time                         \
      --verbose                         \
      --append                          \
      --output="${RESULTS}"             \
      /plink2                           \
          --allow-extra-chr             \
          ${EXTRACT[@]+"${EXTRACT[@]}"} \
          --make-bed                    \
          --new-id-max-allele-len       \
              1000                      \
              missing                   \
          --out                         \
              "${PREFIX}"               \
          --rm-dup                      \
              force-first               \
          --set-all-var-ids             \
              '@:#:$1:$2'               \
          --vcf                         \
              '~{vcf}'


  EXITSTATUS=$?

  (
      exec 1>>"${RESULTS}"

      printf -- '\n\nSTATUS: %d\n\n' "${EXITSTATUS}"

      for FILEPATH in '~{vcf}' "${PREFIX}".*
      do
          stat --format=$'%s\t%n' "${FILEPATH}"
      done
      printf -- '\n\ndf AFTER:\n'
      df --block-size=1 /cromwell_root
      df --human        /cromwell_root
  )
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    disks: "local-disk " + storage + " HDD"
    memory: memory + " GB"
  }

  output {
    File results = results
  }

}
