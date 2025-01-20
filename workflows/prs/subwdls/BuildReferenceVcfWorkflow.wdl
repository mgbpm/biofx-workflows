version 1.0

import "../tasks/HelperTasks.wdl"

workflow BuildReferenceVcf {
  input {
    Array[File]+ weights
    File         pca_variants
    String       workspace
    String       shards_root
    String       work_bucket
    Boolean      nocleanup        = false
    Int          nbatches         = 500
    String       prs_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250122"
  }

  String tmp       = work_bucket + "/.BuildReferenceVcf"
  String workdir   = tmp         + "/work"
  String sentinels = tmp         + "/sentinels"

  if (!nocleanup) {
    scatter (tsv in weights) {
      call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInWeights {
        input:
            tsv        = tsv
          , skipheader = true
      }
    }

    call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInPcaVariants {
      input:
          tsv        = pca_variants
        , skipheader = false
    }
  }

  Array[File] clean_weights      = select_first([RenameChromosomesInWeights.renamed, weights])
  File        clean_pca_variants = select_first([RenameChromosomesInPcaVariants.renamed, pca_variants])

  call HelperTasks.GetTotalSize as FootprintOfWeightsAndPCA {
    input:
        urls         = flatten([clean_weights, [clean_pca_variants]])
      , workspace    = workspace
      , docker_image = prs_docker_image
  }

  call GetRegions {
    input:
        weights       = clean_weights
      , pca_variants  = clean_pca_variants
      , footprint     = FootprintOfWeightsAndPCA.gigabytes
  }

  call HelperTasks.ListShards {
    input:
        source       = shards_root
      , workspace    = workspace
      , docker_image = prs_docker_image
  }

  call PurgeTmp as MaybePurgeTmp {
    input:
        tmp          = tmp
      , workspace    = workspace
      , docker_image = prs_docker_image
  }

  call HelperTasks.MakeBatches {
    input:
        cases    = ListShards.relpaths
      , nbatches = nbatches
  }

  scatter (batch in MakeBatches.batches) {
    call SubsetShards {
      input:
          batch        = batch
        , regions      = GetRegions.regions
        , source       = shards_root
        , target       = workdir
        , sentinels    = sentinels
        , workspace    = workspace
        , docker_image = prs_docker_image
    }
  }

  call HelperTasks.GetTotalSize as FootprintOfSubsettedShards {
    input:
        urls         = [workdir]
      , workspace    = workspace
      , sequencing   = SubsetShards.sequencing
      , docker_image = prs_docker_image
  }

  call ConcatenateShards {
    input:
      basedir      = workdir
    , target       = work_bucket
    , relpaths     = ListShards.relpaths
    , workspace    = workspace
    , storage      = 3 * FootprintOfSubsettedShards.gigabytes + 10
    , docker_image = prs_docker_image
  }

  call PurgeTmp {
    input:
        tmp          = tmp
      , workspace    = workspace
      , sequencing   = ConcatenateShards.sequencing
      , docker_image = prs_docker_image
  }

  output {
    File        reference_vcf        = ConcatenateShards.reference_vcf
    File        reference_tbi        = ConcatenateShards.reference_tbi
    Array[File] renamed_weights      = clean_weights
    File        renamed_pca_variants = clean_pca_variants
    File        regions              = GetRegions.regions
  }
}

# -----------------------------------------------------------------------------

task GetRegions {
  input {
    Array[File]+ weights
    File         pca_variants
    Int          footprint
  }

  String OUTPUTDIR = "OUTPUT"
  String REGIONS   = OUTPUTDIR + "/REGIONS"
  Int    storage   = 20 + 2 * footprint

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  OUTPUTDIR='~{OUTPUTDIR}'
  WEIGHTSFILES=( '~{sep="' '" weights}' )
  PCA='~{pca_variants}'

  mkdir --verbose --parents "${OUTPUTDIR}"

  # ---------------------------------------------------------------------------

  (
    # The perl segment of the following pipeline prints the first line only if
    # it ends in a numeric expression (i.e. the file does not have a headers
    # row)
    for WEIGHTS in "${WEIGHTSFILES[@]}"
    do
        cut --fields=1 "${WEIGHTS}"               \
          | perl -lne '
              BEGIN {
                $::WEIGHT = qr/
                   \b-?
                   (?:\d+|\d*\.\d+|\d+\.\d*)
                   (?:[eE][+-]?\d+)?\b
                /x;
              }
              print if $. > 1 || /$::WEIGHT\s*$/'
    done

    cut --fields=1 "${PCA}"
  )                                               \
    | cut --fields=1,2 --delimiter=:              \
    | sort --unique                               \
    | tr : $'\t'                                  \
    > '~{REGIONS}'
  >>>

  output {
    File regions = REGIONS
  }

  runtime {
    preemptible: 5
    docker     : "ubuntu:21.10"
    disks      : "local-disk ~{storage} HDD"
  }
}

task SubsetShards {
  input {
    Array[String]+ batch
    File           regions
    String         source
    String         target
    String         sentinels
    String         workspace
    String         docker_image
  }

  File firstshard = source + "/" + batch[0]
  Int  storage    = 20 + 3 * ceil(size(firstshard, "GB"))

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  export WORKSPACE='~{workspace}'

  SOURCE="$( mapurl.sh '~{source}' )"
  TARGET="$( mapurl.sh '~{target}' )"
  SENTINELS="$( mapurl.sh '~{sentinels}' )"

  fileexists() {
      local url="${1%/}"
      local parent="$( dirname "${url}" )"
      local basename="$( basename "${url}" )"
      local -a found
      readarray -t found < <(
                              rclone                      \
                                  lsf                     \
                                  --files-only            \
                                  --include="${basename}" \
                                  "${parent}"
                            )

      local nfound="${#found[@]}"
      if (( nfound == 0 ))
      then
          printf -- 'false'
      elif (( nfound == 1 )) && [[ "${found[0]}" == */"${basename}" ]]
      then
          printf -- 'true'
      else
          printf -- 'fileexists: unexpected nfound=%d\n' "${nfound}" >&2
          exit 1
      fi
  }

  subsetshard () {
      local relpath="${1}"
      local sentinel="${SENTINELS}/${relpath}"

      if "$( fileexists "${sentinel}" )"
      then
          return
      fi

      local shard="${SOURCE}/${relpath}"

      local inputshard="$( mktemp )"
      local outputshard="$( mktemp )"

      rclone copyto "${SOURCE}/${relpath}"     "${inputshard}"
      rclone copyto "${SOURCE}/${relpath}.tbi" "${inputshard}.tbi"

      touch --date '1 Jan 1970 00:00:00 -0000' "${inputshard}"
      touch --date '1 Jan 1970 00:00:01 -0000' "${inputshard}.tbi"

      bcftools                          \
          view                          \
          --no-version                  \
          --output-type v               \
          --regions-file '~{regions}'   \
          "${inputshard}"               \
        | bcftools                      \
              norm                      \
              --multiallelics -any      \
              --no-version              \
              --output-type v           \
        | bcftools                      \
              annotate                  \
              --remove 'INFO,FORMAT'    \
              --no-version              \
              --output-type v           \
        | bcftools                      \
              view                      \
              --no-version              \
              --output-type z           \
              --output "${outputshard}"

      bcftools             \
          index            \
          --force          \
          --tbi            \
          "${outputshard}"

      rclone copyto "${outputshard}"     "${TARGET}/${relpath}"
      rclone copyto "${outputshard}.tbi" "${TARGET}/${relpath}.tbi"

      rclone touch "${sentinel}"

      rm --force "${inputshard}"* "${outputshard}"*
  }

  BATCH=( '~{sep="' '" batch}' )

  for RELPATH in "${BATCH[@]}"
  do
      subsetshard "${RELPATH}"
  done
  >>>

  output {
    Boolean sequencing = true
  }

  runtime {
    preemptible: 5
    docker     : docker_image
    disks      : "local-disk ~{storage} HDD"
  }
}

task ConcatenateShards {
  input {
    String basedir
    String target
    File   relpaths
    String workspace
    Int    storage
    String docker_image

    Int    memory       = 8
  }

  String OUTPUTDIR = "OUTPUT"
  String REFERENCE = OUTPUTDIR + "/reference.vcf.gz"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  set -o xtrace

  # export RCLONE_LOG_LEVEL=DEBUG
  # export RCLONE_STATS_LOG_LEVEL=DEBUG

  # typeset -p >&2
  # rclone version >&2

  free --bytes >&2
  free --human >&2

  export WORKSPACE='~{workspace}'

  BASEDIR="$( mapurl.sh '~{basedir}' )"
  TARGET="$( mapurl.sh '~{target}' )"

  LOCALBASEDIR="$( mktemp --directory )"

  rclone copy "${BASEDIR}" "${LOCALBASEDIR}"

  find "${LOCALBASEDIR}" -type f

  SHARDS="$( mktemp )"
  perl -lpe "s@^@${LOCALBASEDIR}/@" '~{relpaths}' > "${SHARDS}"

  mkdir --verbose --parents '~{OUTPUTDIR}'

  bcftools                    \
      concat                  \
      --file-list="${SHARDS}" \
      --no-version            \
      --naive                 \
    > '~{REFERENCE}'

  rm --recursive --force "${LOCALBASEDIR}"

  bcftools             \
      index            \
      --force          \
      --tbi            \
      '~{REFERENCE}'

  BASENAME="$( basename '~{REFERENCE}' )"

  for EXTENSION in '' '.tbi'
  do
      rclone copyto "~{REFERENCE}${EXTENSION}" "${TARGET}/${BASENAME}${EXTENSION}"
  done

  >>>

  output {
    File    reference_vcf = REFERENCE
    File    reference_tbi = REFERENCE + ".tbi"
    Boolean sequencing    = true
  }

  runtime {
    preemptible : 1
    disks       : "local-disk " + storage + " HDD"
    memory      : memory + " GB"
    docker      : docker_image
  }
}

task PurgeTmp {
  input {
    String tmp
    String workspace
    String docker_image
    String sequencing   = "ignored"
  }

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  # set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  export WORKSPACE='~{workspace}'
  TMP="$( mapurl.sh '~{tmp}' )"
  rclone purge "${TMP}"
  >>>
  runtime {
    preemptible: 5
    docker     : docker_image
  }
}
