version 1.0

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/workflows/prs/tasks/HelperTasks.wdl" as HelperTasks
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/workflows/prs/tasks/ScoringTasks.wdl" as ScoringTasks

workflow PrsInputPrep {
  input {
    Array[File] weights_files
    File        pca_variants
    String      workspace
    String      source
    String      target
    Int         nbatches         = 500
    Boolean     resuming         = false
    Boolean     norename         = false
    Array[File] query_vcfs
    String      prs_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515"
  }

  String tmp              = target + "/.PreparePrsInputs"
  String workdir          = tmp    + "/work"
  String sentinels        = tmp    + "/sentinels"

  if (! norename) {
    scatter (weights in weights_files) {
      call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInWeights {
        input:
            tsv        = weights
          , skipheader = true
      }
    }

    call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInPcaVariants {
      input:
          tsv        = pca_variants
        , skipheader = false
    }
  }

  Array[File] weights_files_ = select_first([RenameChromosomesInWeights.renamed,
                                             weights_files])
  File        pca_variants_  = select_first([RenameChromosomesInPcaVariants.renamed,
                                             pca_variants])

  call GetTotalSize as FootprintOfWeightsAndPCA {
    input:
        urls         = flatten([weights_files_, [pca_variants_]])
      , workspace    = workspace
      , docker_image = prs_docker_image
  }

  call GetRegions {
    input:
        weights_files = weights_files_
      , pca_variants  = pca_variants_
      , footprint     = FootprintOfWeightsAndPCA.gigabytes
  }

  call HelperTasks.ListShards {
    input:
        source       = source
      , workspace    = workspace
      , docker_image = prs_docker_image
  }

  if (resuming) {
    call FetchSentinels {
      input:
          basedir      = sentinels
        , workspace    = workspace
        , docker_image = prs_docker_image
    }
  }

  if (!resuming) {
    call PurgeTmp as MaybePurgeTmp {
      input:
          tmp          = tmp
        , workspace    = workspace
        , docker_image = prs_docker_image
    }
  }

  call HelperTasks.MakeBatches {
    input:
        cases    = ListShards.relpaths
      , nbatches = nbatches
      , exclude  = select_first([FetchSentinels.sentinels, []])
  }

  scatter (batch in MakeBatches.batches) {
    call SubsetShards {
      input:
          batch        = batch
        , regions      = GetRegions.regions
        , source       = source
        , target       = workdir
        , sentinels    = sentinels
        , workspace    = workspace
        , docker_image = prs_docker_image
    }
  }

  call GetTotalSize as FootprintOfSubsettedShards {
    input:
        urls         = [workdir]
      , workspace    = workspace
      , sequencing   = SubsetShards.sequencing
      , docker_image = prs_docker_image
  }

  call ConcatenateShards {
    input:
      basedir      = workdir
    , target       = target
    , relpaths     = ListShards.relpaths
    , workspace    = workspace
    , storage      = 3 * FootprintOfSubsettedShards.gigabytes + 10
    , docker_image = prs_docker_image
  }

  call HelperTasks.GetBaseMemory as GetMemoryForReference {
    input:
        vcf = ConcatenateShards.reference_vcf
  }

  call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
    input:
        vcf = ConcatenateShards.reference_vcf
      , mem = GetMemoryForReference.gigabytes
  }

  scatter (query_vcf in query_vcfs) {
    call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
      input:
          vcf = query_vcf
    }

    if (! norename) {
      call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
        input:
            vcf = query_vcf
      }
    }

    File query_vcf_ = select_first([RenameChromosomesInQueryVcf.renamed,
                                    query_vcf])

    call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
      input:
          vcf = query_vcf_
        , mem = GetMemoryForQueryFromVcf.gigabytes
    }
  }

  call GetTotalSize as FootprintOfVariantFiles {
    input:
        urls         = ExtractQueryVariants.ids
      , workspace    = workspace
      , docker_image = prs_docker_image
  }

  call HelperTasks.Union {
    input:
        lists   = ExtractQueryVariants.ids
      , storage = FootprintOfVariantFiles.gigabytes
  }

  call HelperTasks.Intersection {
    input:
        lists   = ExtractQueryVariants.ids
      , storage = FootprintOfVariantFiles.gigabytes
  }

  call MaybeTrimPcaVariants {
    input:
        pca_variants = pca_variants_
      , reference    = ExtractReferenceVariants.ids
      , union        = Union.result
      , intersection = Intersection.result
  }

  call PurgeTmp {
    input:
        tmp          = tmp
      , workspace    = workspace
      , sequencing   = MaybeTrimPcaVariants.sequencing
      , docker_image = prs_docker_image
  }

  output {
    File        regions               = GetRegions.regions
    File        kept_pca_variants     = select_first([MaybeTrimPcaVariants.kept_pca_variants,
                                                      pca_variants_])
    Array[File] renamed_weights_files = weights_files_
    Array[File] renamed_query_vcfs    = query_vcf_
    File        reference_vcf         = ConcatenateShards.reference_vcf
    File        reference_tbi         = ConcatenateShards.reference_tbi
  }
}

# -----------------------------------------------------------------------------

task GetTotalSize {
  input {
    Array[String]+  urls
    String          workspace
    String          docker_image
    Array[Boolean]+ sequencing = [false]  # ignored!
  }

  String OUTPUTDIR = "OUTPUT"
  String GIGABYTES = OUTPUTDIR + "/GIGABYTES"

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  export WORKSPACE='~{workspace}'
  OUTPUTDIR='~{OUTPUTDIR}'

  mkdir --verbose --parents "${OUTPUTDIR}"

  TOTALBYTES=0
  for URL in '~{sep="' '" urls}'
  do
      BYTES="$(
                rclone size "$( mapurl.sh "${URL}" )"  \
                  | tail --lines=1                     \
                  | perl -lpe 's/^.*?(\d+) Byte.*/$1/'
              )"
      TOTALBYTES=$(( TOTALBYTES + BYTES ))
  done

  perl -e "printf qq(%f\n), ${TOTALBYTES}/2**30" > '~{GIGABYTES}'

  printf -- '\nTotal size: %.1f GiB\n' "$( cat '~{GIGABYTES}' )"
  >>>

  output {
    Int gigabytes = ceil(read_float(GIGABYTES))
  }

  runtime {
    preemptible: 5
    docker     : docker_image
  }
}

task GetRegions {
  input {
    Array[File]+ weights_files
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
  WEIGHTSFILES=( '~{sep="' '" weights_files}' )
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

task FetchSentinels {
  input {
    String basedir
    String workspace
    String docker_image
  }

  String SENTINELS = "output/SENTINELS"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  typeset -p >&2
  rclone version >&2

  mkdir --verbose --parents "$( dirname '~{SENTINELS}' )"

  export WORKSPACE='~{workspace}'
  BASEDIR="$( mapurl.sh '~{basedir}' )"

  rclone             \
      lsf            \
      --recursive    \
      --files-only   \
      "${BASEDIR}"   \
    > '~{SENTINELS}'

  >>>

  output {
    Array[String] sentinels = read_lines(SENTINELS)
  }

  runtime {
    preemptible: 5
    docker     : docker_image
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

task MaybeTrimPcaVariants {
  input {
    File    pca_variants
    File    reference
    File    union
    File    intersection
    Boolean debug        = false
  }

  String        WORKDIR                 = "WORK"
  String        ARCHIVE                 = WORKDIR + ".tgz"
  String        OUTPUTDIR               = "OUTPUT"
  String        KEPT_PCA_VARIANTS       = OUTPUTDIR + "/kept_pca_variants.tsv"
  String        WARNINGS                = OUTPUTDIR + "/WARNINGS"
  Array[String] WORKFILES               = [
                                              WORKDIR + "/PCA"
                                            , WORKDIR + "/UNION"
                                            , WORKDIR + "/INTERSECTION"
                                            , WORKDIR + "/REFERENCE"
                                            , WORKDIR + "/PU"
                                            , WORKDIR + "/PI"
                                            , WORKDIR + "/PR"
                                            , WORKDIR + "/NIXQ"
                                            , WORKDIR + "/NIXR"
                                            , WORKDIR + "/NIX"
                                            , WORKDIR + "/TEMP_PCA"
                                            , WORKDIR + "/WANTED_PCA"
                                          ]
  Int           multiplier          = if debug then 20 else 10
  Int           storage             = 30 + multiplier * ceil(  size(pca_variants, "GB")
                                                             + size(union       , "GB")
                                                             + size(intersection, "GB")
                                                             + size(reference   , "GB"))

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  alert() {
      local tag="${1}"
      local message="${2}"
      printf -- '%s: %s\n' "${tag}" "${message}"
  }

  error() {
      local message="${1}"
      alert 'ERROR' "${message}" >&2
      exit 1
  }

  warning() {
      local message="${1}"
      alert 'WARNING' "${message}" >> '~{WARNINGS}'
  }

  # ---------------------------------------------------------------------------

  printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
  printf -- 'INITIAL STORAGE UTILIZATION:\n'
  df --human
  printf -- '\n'

  # ---------------------------------------------------------------------------

  WORKDIR='~{WORKDIR}'
  OUTPUTDIR='~{OUTPUTDIR}'
  UNSORTEDPCA='~{pca_variants}'    # pca_variants file (no header row)
  UNSORTEDREFERENCE='~{reference}' # CHROM:POS:REF:ALT from reference vcf

  # ---------------------------------------------------------------------------

  mkdir --verbose --parents "${WORKDIR}"
  mkdir --verbose --parents "${OUTPUTDIR}"

  WIDTH=$(( ${#WORKDIR} + 10 ))
  TEMPLATE="Creating %-${WIDTH}s ... "
  unset WIDTH

  PCA="${WORKDIR}/PCA"
  REFERENCE="${WORKDIR}/REFERENCE"
  UNION="${WORKDIR}/UNION"
  INTERSECTION="${WORKDIR}/INTERSECTION"

  printf -- "${TEMPLATE}" "${PCA}"
  sort --unique "${UNSORTEDPCA}" > "${PCA}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${REFERENCE}"
  sort --unique "${UNSORTEDREFERENCE}" > "${REFERENCE}"
  printf -- 'done\n'

  # ---------------------------------------------------------------------------

  ln --symbolic --verbose '~{union}'        "${UNION}"
  ln --symbolic --verbose '~{intersection}' "${INTERSECTION}"

  PU="${WORKDIR}/PU"
  printf -- "${TEMPLATE}" "${PU}"
  comm -1 -2 "${UNION}" "${PCA}" > "${PU}"
  printf -- 'done\n'

  if [[ ! -s ${PU} ]]
  then
      error 'No variant in the PCA variants file is mentioned in the query VCFs.'
  fi

  PI="${WORKDIR}/PI"
  printf -- "${TEMPLATE}" "${PI}"
  comm -1 -2 "${INTERSECTION}" "${PCA}" > "${PI}"
  printf -- 'done\n'

  if [[ ! -s ${PI} ]]
  then
      warning 'No variant in the PCA variants file is mentioned in all the query VCFs.'
  fi

  PR="${WORKDIR}/PR"
  printf -- "${TEMPLATE}" "${PR}"
  comm -1 -2 "${REFERENCE}" "${PCA}" > "${PR}"
  printf -- 'done\n'

  if [[ ! -s ${PR} ]]
  then
      error 'No variant in the PCA variants file is mentioned in the reference VCF.'
  fi

  NPR=$( wc --lines < "${PR}" )
  NPCS=20
  if (( NPR < 2 * NPCS + 1 ))
  then

      # NB: The flashpca program, invoked by PCATasks.PerformPCA, will fail if
      #
      #     min(NVARIANTS, NSAMPLES) < 2 * NPCS + 1
      #
      # ...where NVARIANTS and NSAMPLES are the numbers of variants (rows) and
      # samples (data columns) in the reference VCF, and NPCS is the value of
      # the command's --ndim flag (currently hard-coded as 20; see
      # PCATasks.PerformPCA).  (The inequality above is derived from line 623 of
      # https://github.com/gabraham/flashpca/blob/b8044f13607a072125828547684fde8b081d6191/flashpca.cpp .)
      # In particular, flashpca will fail when NPR = NVARIANTS is less than
      # 2 * NPCS + 1.  Hence the error condition detected here.

      error 'The reference VCF mentions too few PCA variants for the PCA calculation.'
  fi

  NIXQ="${WORKDIR}/NIXQ"
  printf -- "${TEMPLATE}" "${NIXQ}"
  comm -2 -3 "${PU}" "${PR}" > "${NIXQ}"
  printf -- 'done\n'

  NIXR="${WORKDIR}/NIXR"
  printf -- "${TEMPLATE}" "${NIXR}"
  comm -2 -3 "${PR}" "${PI}" > "${NIXR}"
  # Equivalently, one could do this:
  # comm -2 -3 "${PR}" "${INTERSECTION}" > "${NIXR}"
  printf -- 'done\n'

  NIX="${WORKDIR}/NIX"
  printf -- "${TEMPLATE}" "${NIX}"
  sort "${NIXQ}" "${NIXR}" > "${NIX}"
  printf -- 'done\n'

  if [[ ! -s ${NIX} ]]
  then
      exit 0
  fi

  (
      sortvariants() {
          local variants="${1}"
          sort                    \
              --field-separator=: \
              --key=1,1V          \
              --key=2,2n          \
              --key=3,3           \
              --key=4,4           \
              "${variants}"
      }

      SEPARATOR=''
      maybewarn() {
          local nixx="${1}"
          local label0="${2}"
          local label1
          if [[ ${label0} == query ]]
          then
              label1=reference
          else
              label1=query
          fi
          if [[ -s ${nixx} ]]
          then
              printf -- "${SEPARATOR}"
              if [[ -z "${SEPARATOR}" ]]
              then
                  SEPARATOR=$'\n'
              fi
              printf -- 'WARNING: '
              printf -- 'the following PCA variants from the '
              printf -- '%s VCF do not appear in the %s VCF, '  \
                        "${label0}" "${label1}"
              printf -- 'and will be removed from the '
              printf -- 'analysis.\n'
              sortvariants "${nixx}"
          fi
      }

      exec >>'~{WARNINGS}'
      maybewarn "${NIXQ}" query
      maybewarn "${NIXR}" reference
  )

  # ---------------------------------------------------------------------------

  TEMP_PCA="${WORKDIR}/TEMP_PCA"
  printf -- "${TEMPLATE}" "${TEMP_PCA}"
  comm -2 -3 "${PCA}" "${NIX}" > "${TEMP_PCA}"
  printf -- 'done\n'

  WANTED_PCA="${WORKDIR}/WANTED_PCA"
  printf -- "${TEMPLATE}" "${WANTED_PCA}"
  join                                                \
      -t $'\t'                                        \
      "${TEMP_PCA}"                                   \
      <(
         perl -lpe '$_ = qq($_\t$.)' "${UNSORTEDPCA}" \
           | sort                                     \
                 --field-separator=$'\t'              \
                 --key=1,1
       )                                              \
    | sort                                            \
          --field-separator=$'\t'                     \
          --key=2,2n                                  \
    | cut --fields=1                                  \
    > "${WANTED_PCA}"
  printf -- 'done\n'

  mkdir --verbose --parents "$( dirname '~{KEPT_PCA_VARIANTS}' )"
  cp --verbose "${WANTED_PCA}" '~{KEPT_PCA_VARIANTS}'

  # ---------------------------------------------------------------------------

  printf -- '\n\n## WORKDIR:\n'
  find -L '~{WORKDIR}' -type f | xargs ls -ltr

  (
      exec >&2
      printf -- '\n\n### LINE COUNTS:\n'
      find -L '~{WORKDIR}' -type f | xargs ls -1tr | xargs wc --lines
  )

  # ---------------------------------------------------------------------------

  if ~{if debug then "true" else "false"}
  then
      tar -cvzf '~{ARCHIVE}' '~{WORKDIR}'
  fi

  # ---------------------------------------------------------------------------

  printf -- 'FINAL STORAGE UTILIZATION:\n'
  df --human

  # ---------------------------------------------------------------------------

  if ~{if !debug then "true" else "false"}
  then
      rm -rf '~{WORKDIR}'
  fi

  # ---------------------------------------------------------------------------
  >>>

  output {
    File?        kept_pca_variants    = KEPT_PCA_VARIANTS
    File?        warnings             = WARNINGS

    Array[File?] workfiles            = WORKFILES
    File?        workfiles_tgz        = ARCHIVE
    Boolean      sequencing           = true
  }

  runtime {
    preemptible: 5
    disks      : "local-disk ~{storage} HDD"
    docker     : "ubuntu:21.10"
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
