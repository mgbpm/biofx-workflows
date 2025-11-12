version 1.0

import "../tasks/HelperTasks.wdl"
import "../tasks/ScoringTasks.wdl"

workflow PreparePrsInputs {
  input {
    Array[File] variant_weights
    File        pca_variants
    String      workspace
    String      source
    String      target
    Int         nbatches            = 500
    Boolean     resuming            = false
    Boolean     norename            = false
    Array[File] query_vcfs
    String      prs_docker_image    = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515"
    String      ubuntu_docker_image = "ubuntu:latest"
  }

  String tmp              = target + "/.PreparePrsInputs"
  String workdir          = tmp    + "/work"
  String sentinels        = tmp    + "/sentinels"

  if (!norename) {
    scatter (weights in variant_weights) {
      call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInWeights {
        input:
          tsv = weights,
          skipheader = true
      }
    }

    call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInPcaVariants {
      input:
        tsv = pca_variants,
        skipheader = false
    }
  }

  Array[File] variant_weights_ = select_first([RenameChromosomesInWeights.renamed, variant_weights])
  File pca_variants_  = select_first([RenameChromosomesInPcaVariants.renamed, pca_variants])

  call GetTotalSize as FootprintOfWeightsAndPCA {
    input:
      urls = flatten([variant_weights_, [pca_variants_]]),
      workspace = workspace,
      docker_image = prs_docker_image
  }

  call GetRegions {
    input:
      variant_weights = variant_weights_,
      pca_variants = pca_variants_,
      footprint = FootprintOfWeightsAndPCA.gigabytes
  }

  call HelperTasks.ListShards {
    input:
      source = source,
      workspace = workspace,
      docker_image = prs_docker_image
  }

  if (resuming) {
    call FetchSentinels {
      input:
        basedir = sentinels,
        workspace = workspace,
        docker_image = prs_docker_image
    }
  }

  if (!resuming) {
    call PurgeTmp as MaybePurgeTmp {
      input:
        tmp = tmp,
        workspace = workspace,
        docker_image = prs_docker_image
    }
  }

  call HelperTasks.MakeBatches {
    input:
      cases = ListShards.relpaths,
      nbatches = nbatches,
      exclude  = select_first([FetchSentinels.sentinels, []])
  }

  scatter (batch in MakeBatches.batches) {
    call SubsetShards {
      input:
        batch = batch,
        regions = GetRegions.regions,
        source = source,
        target = workdir,
        sentinels = sentinels,
        workspace = workspace,
        docker_image = prs_docker_image
    }
  }

  call GetTotalSize as FootprintOfSubsettedShards {
    input:
      urls = [workdir],
      workspace = workspace,
      sequencing = SubsetShards.sequencing,
      docker_image = prs_docker_image
  }

  call ConcatenateShards {
    input:
      basedir = workdir,
      target = target,
      relpaths = ListShards.relpaths,
      workspace = workspace,
      storage = 3 * FootprintOfSubsettedShards.gigabytes + 10,
      docker_image = prs_docker_image
  }

  call HelperTasks.GetBaseMemory as GetMemoryForReference {
    input:
        vcf = ConcatenateShards.reference_vcf
  }

  call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
    input:
      vcf = ConcatenateShards.reference_vcf,
      mem_size = GetMemoryForReference.gigabytes
  }

  scatter (query_vcf in query_vcfs) {
    call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
      input:
          vcf = query_vcf
    }

    if (!norename) {
      call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
        input:
            vcf = query_vcf
      }
    }

    File query_vcf_ = select_first([RenameChromosomesInQueryVcf.renamed, query_vcf])

    call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
      input:
        vcf = query_vcf_,
        mem_size = GetMemoryForQueryFromVcf.gigabytes
    }
  }

  call GetTotalSize as FootprintOfVariantFiles {
    input:
      urls = ExtractQueryVariants.ids,
      workspace = workspace,
      docker_image = prs_docker_image
  }

  call HelperTasks.Union {
    input:
      lists = ExtractQueryVariants.ids,
      addldisk = FootprintOfVariantFiles.gigabytes
  }

  call HelperTasks.Intersection {
    input:
      lists = ExtractQueryVariants.ids,
      addldisk = FootprintOfVariantFiles.gigabytes
  }

  call HelperTasks.TrimPcaToSubset as TrimVariants {
    input:
        pca_variants = pca_variants_
      , reference    = ExtractReferenceVariants.ids
      , union        = Union.result
      , intersection = Intersection.result
  }

  output {
    File        regions                 = GetRegions.regions
    File?       kept_pca_variants       = TrimVariants.kept_pca_variants
    Array[File] renamed_variant_weights = variant_weights_
    Array[File] renamed_query_vcfs      = query_vcf_
    File        reference_vcf           = ConcatenateShards.reference_vcf
    File        reference_tbi           = ConcatenateShards.reference_tbi
  }
}

task GetTotalSize {
  input {
    Array[String]+  urls
    String          workspace
    String          docker_image
    Array[Boolean]+ sequencing  = [false]  # ignored!
    Int             preemptible = 5
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
    docker: "~{docker_image}"
    preemptible: preemptible
  }
}

task GetRegions {
  input {
    Array[File]+ variant_weights
    File         pca_variants
    String       docker_image = "ubuntu:21.10"
    Int          footprint
    Int          preemptible  = 5
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
    WEIGHTSFILES=( '~{sep="' '" variant_weights}' )
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
    docker: "~{docker_image}"
    disks: "local-disk ~{storage} HDD"
    preemptible: preemptible
  }
}

task FetchSentinels {
  input {
    String basedir
    String workspace
    String docker_image
    Int    preemptible = 5
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
    docker: "~{docker_image}"
    preemptible: preemptible
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
    Int            preemptible = 5
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
    docker : "~{docker_image}"
    disks : "local-disk ~{storage} HDD"
    preemptible: preemptible
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
    Int    mem_size    = 8
    Int    preemptible = 1
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
    docker: "~{docker_image}"
    disks: "local-disk ~{storage} HDD"
    memory: "~{mem_size} GB"
    preemptible : preemptible
  }
}

task PurgeTmp {
  input {
    String tmp
    String workspace
    String docker_image
    String sequencing  = "ignored"
    Int    preemptible = 5
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
    docker: "~{docker_image}"
    preemptible: preemptible
  }
}
