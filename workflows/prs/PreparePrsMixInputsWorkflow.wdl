version 1.0

import "ScoringTasks.wdl"
import "HelperTasks.wdl"

workflow PreparePrsMixInputs {
  input {
    Array[File] weights_files
    File        pca_variants
    String      workspace
    String      source
    String      target
    Int         nbatches  = 500
    Boolean     resuming  = false
  }

  String tmp       = target + "/.PreparePrsInputs"
  String workdir   = tmp    + "/work"
  String sentinels = tmp    + "/sentinels"

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

  call GetTotalSize as FootprintOfWeightsAndPCA {
    input:
        urls      = flatten([RenameChromosomesInWeights.renamed,
                             [RenameChromosomesInPcaVariants.renamed]])
      , workspace = workspace
  }

  call GetReferenceRegions {
    input:
        weights_files = RenameChromosomesInWeights.renamed
      , pca_variants  = RenameChromosomesInPcaVariants.renamed
      , footprint     = FootprintOfWeightsAndPCA.gigabytes
  }

  # FIXME: the image for this should a minimal image + rclone
  call HelperTasks.ListShards {
    input:
        source       = source
      , workspace    = workspace
      , docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/sharding:0.0.1"
  }

  if (resuming) {
    call FetchSentinels {
      input:
          basedir   = sentinels
        , workspace = workspace
    }
  }

  if (!resuming) {
    call PurgeTmp as MaybePurgeTmp {
      input:
          tmp       = tmp
        , workspace = workspace
    }
  }

  call HelperTasks.MakeBatches {
    input:
        cases    = ListShards.relpaths
      , nbatches = nbatches
      , exclude  = select_first([FetchSentinels.sentinels, []])
  }

  # scatter (batch in MakeBatches.batches) {
  #   call SubsetShards {
  #     input:
  #         batch     = batch
  #       , source    = source
  #       , target    = workdir
  #       , regions   = regions
  #       , workspace = workspace
  #   }
  # }
  #
  # call GetTotalSize as FootprintOfSubsettedShards {
  #   input:
  #       urls       = [workdir]
  #     , sequencing = SubsetShards.sequencing
  # }
  #
  # call ConcatenateShards {
  #   input:
  #     docker_image = docker_image
  #   , storage      = 3 * FootprintOfSubsettedShards.footprint + 10
  #   , basedir      = basedir
  #   , shards       = spec.shards
  #   , target       = target
  # }
  #
  # call MakeReferenceVcf {
  #   input:
  #       regions    = GetReferenceRegions.regions
  #     , shardsbase = shardsbase
  # }
  #
  # call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
  #   input:
  #       vcf = MakeReferenceVcf
  #     , mem = GetMemoryForReference.gigabytes
  # }
  #
  # scatter (vcf in queryvcfs) {
  #   call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
  #     input:
  #         vcf = query_file
  #   }
  #
  #   call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
  #     input:
  #         vcf = query_file
  #   }
  #
  #   call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
  #     input:
  #         vcf = RenameChromosomesInQueryVcf.renamed
  #       , mem = GetMemoryForQueryFromVcf.gigabytes
  #   }
  # }
  #
  # call GetTotalSize as FootprintOfVariantFiles {
  #   input:
  #       urls = ExtractQueryVariants.ids
  # }
  #
  # call HelperTasks.Union {
  #   input:
  #       lists   = ExtractQueryVariants.ids
  #     , storage = FootprintOfVariantFiles.gigabytes
  # }
  #
  # call HelperTasks.Intersection {
  #   input:
  #       lists   = ExtractQueryVariants.ids
  #     , storage = FootprintOfVariantFiles.gigabytes
  # }
  #
  # call MaybeTrimPcaVariants {
  #   input:
  #       pca_variants = RenameChromosomesInPcaVariants.renamed
  #     , reference    = ExtractReferenceVariants.ids
  #     , union        = Union.result
  #     , intersection = Intersection.result
  # }
  #
  # call PurgeTmp {
  #   input:
  #       tmp              = tmp
  #     , sequencing_dummy = MaybeTrimPcaVariants.sequencing_dummy
  # }

}

# -----------------------------------------------------------------------------

task GetTotalSize {
  input {
    Array[String]+  urls
    String          workspace
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
    # FIXME: the image for this should a minimal image + rclone
    docker     : "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/sharding:0.0.1"
  }
}

task GetReferenceRegions {
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
    docker: "ubuntu:21.10"
    disks : "local-disk ~{storage} HDD"
  }
}

task FetchSentinels {
  input {
    String basedir
    String workspace
  }

  command <<<
  >>>

  output {
    Array[String] sentinels = []
  }

  runtime {
    # FIXME: the image for this should a minimal image + rclone
    docker     : "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/sharding:0.0.1"
  }
}

# task FetchSentinels {
#   input {
#     String basedir
#     String workspace
#   }
#
#   String SENTINELS = "output/SENTINELS"
#
#   command <<<
#   set -o errexit
#   set -o pipefail
#   set -o nounset
#   # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
#   # set -o xtrace
#
#   typeset -p >&2
#   rclone version >&2
#
#   mkdir --verbose --parents "$( dirname '~{SENTINELS}' )"
#
#   export WORKSPACE='~{workspace}'
#   BASEDIR="$( mapurl.sh '~{basedir}' )"
#
#   rclone             \
#       lsf            \
#       --recursive    \
#       --files-only   \
#       "${BASEDIR}"   \
#     > '~{SENTINELS}'
#
#   >>>
#
#   output {
#     Array[String] sentinels = read_lines(SENTINELS)
#   }
#
#   runtime {
#     # FIXME: the image for this should a minimal image + rclone
#     docker     : "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/sharding:0.0.1"
#     preemptible: 5
#   }
# }

# task SubsetShards {
#   input {
#     Array[String]+ batch
#     String         source
#     String         target
#     File           regions
#     String         workspace
#   }
#
#   File firstshard = source + "/" + batch[0]
#   Int  storage   = 20 + 2 * ceil(size(firstshard))
#
#   command <<<
#   set -o errexit
#   set -o pipefail
#   set -o nounset
#   set -o xtrace
#   # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
#
#   export WORKSPACE='~{workspace}'
#
#   subset_shard () {
#       local relpath="${1}"
#       local source="${2}"
#       local target="${3}"
#       local shard="${source}/${relpath}"
#
#       local inputshard="$( mktemp )"
#       local outputshard="$( mktemp )"
#
#       rclone copyto "${source}/${relpath}"     "${inputshard}"
#       rclone copyto "${source}/${relpath}.tbi" "${inputshard}.tbi"
#
#       touch --date '1 Jan 1970 00:00:00 -0000' "${inputshard}"
#       touch --date '1 Jan 1970 00:00:01 -0000' "${inputshard}.tbi"
#
#       bcftools                          \
#           view                          \
#           --no-version                  \
#           --output-type v               \
#           --regions-file '~{regions}'   \
#           "${inputshard}"               \
#         | bcftools                      \
#               norm                      \
#               --multiallelics -any      \
#               --no-version              \
#               --output-type v           \
#         | bcftools                      \
#               annotate                  \
#               --remove 'INFO,FORMAT'    \
#               --no-version              \
#               --output-type v           \
#         | bcftools                      \
#               view                      \
#               --no-version              \
#               --output-type z           \
#               --output "${outputshard}"
#
#       bcftools             \
#           index            \
#           --force          \
#           --tbi            \
#           "${outputshard}"
#
#       rclone copyto "${outputshard}"     "${target}/${relpath}"
#       rclone copyto "${outputshard}.tbi" "${target}/${relpath}.tbi"
#
#       rm --force "${inputshard}"* "${outputshard}"*
#   }
#
#   SOURCE="$( mapurl.sh '~{source}' )"
#   TARGET="$( mapurl.sh '~{target}' )"
#   BATCH=( '~{sep="' '" batch}' )
#
#   for RELPATH in "${BATCH[@]}"
#   do
#       subset_shard "${RELPATH}" "${SOURCE}" "${TARGET}"
#   done
#   >>>
#
#   output {
#     Boolean sequencing = true
#   }
#
#   runtime {
#     # FIXME: the image for this should a minimal image + rclone
#     docker: "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/sharding:0.0.1"
#     disks : "local-disk ~{storage} HDD"
#   }
# }
#
# task MakeReferenceVcf {
#   input {
#     File   weights
#     File   pca_variants
#     String shardsbase
#   }
# }

# task MaybeTrimPcaVariants {
#   input {
#     File    pca_variants
#     File    reference
#     File    union
#     File    intersection
#     Boolean debug        = false
#   }
#
#   String        WORKDIR                 = "WORK"
#   String        ARCHIVE                 = WORKDIR + ".tgz"
#   String        OUTPUTDIR               = "OUTPUT"
#   String        KEPT_PCA_VARIANTS       = OUTPUTDIR + "/kept_pca_variants.txt"
#   String        WARNINGS                = OUTPUTDIR + "/WARNINGS"
#   Array[String] WORKFILES               = [
#                                               WORKDIR + "/PCA"
#                                             , WORKDIR + "/QUERY"
#                                             , WORKDIR + "/UNION"
#                                             , WORKDIR + "/INTERSECTION"
#                                             , WORKDIR + "/REFERENCE"
#                                             , WORKDIR + "/PQ"
#                                             , WORKDIR + "/PR"
#                                             , WORKDIR + "/NIXQ"
#                                             , WORKDIR + "/NIXR"
#                                             , WORKDIR + "/NIX"
#                                             , WORKDIR + "/TEMP_PCA"
#                                             , WORKDIR + "/WANTED_PCA"
#                                           ]
#   Int           multiplier          = if debug then 20 else 10
#   Int           storage             = 30 + multiplier * ceil(  size(pca_variants, "GB")
#                                                              + size(query       , "GB")
#                                                              + size(reference   , "GB"))
#
#   command <<<
#   set -o pipefail
#   set -o errexit
#   set -o nounset
#   # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
#   # set -o xtrace
#
#   # ---------------------------------------------------------------------------
#
#   error() {
#       local message="${1}"
#       printf -- '%s\n' "${message}" >&2
#       exit 1
#   }
#
#   # ---------------------------------------------------------------------------
#
#   printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
#   printf -- 'INITIAL STORAGE UTILIZATION:\n'
#   df --human
#   printf -- '\n'
#
#   # ---------------------------------------------------------------------------
#
#   WORKDIR='~{WORKDIR}'
#   OUTPUTDIR='~{OUTPUTDIR}'
#   UNSORTEDPCA='~{pca_variants}'    # pca_variants file (no header row)
#   UNSORTEDREFERENCE='~{reference}' # CHROM:POS:REF:ALT from reference vcf
#   # UNSORTEDQUERY='~{query}'         # CHROM:POS:REF:ALT from query vcf
#
#   # ---------------------------------------------------------------------------
#
#   mkdir --verbose --parents "${WORKDIR}"
#   mkdir --verbose --parents "${OUTPUTDIR}"
#
#   WIDTH=$(( ${#WORKDIR} + 10 ))
#   TEMPLATE="Creating %-${WIDTH}s ... "
#   unset WIDTH
#
#   PCA="${WORKDIR}/PCA"
#   REFERENCE="${WORKDIR}/REFERENCE"
#   # QUERY="${WORKDIR}/QUERY"
#   # UNION="${WORKDIR}/UNION"
#   # INTERSECTION="${WORKDIR}/INTERSECTION"
#
#   printf -- "${TEMPLATE}" "${PCA}"
#   sort --unique "${UNSORTEDPCA}" > "${PCA}"
#   printf -- 'done\n'
#
#   # printf -- "${TEMPLATE}" "${QUERY}"
#   # sort --unique "${UNSORTEDQUERY}" > "${QUERY}"
#   # printf -- 'done\n'
#
#   printf -- "${TEMPLATE}" "${REFERENCE}"
#   sort --unique "${UNSORTEDREFERENCE}" > "${REFERENCE}"
#   printf -- 'done\n'
#
#   # ---------------------------------------------------------------------------
#
#   PQ="${WORKDIR}/PQ"
#   printf -- "${TEMPLATE}" "${PQ}"
#   comm -1 -2 "${QUERY}" "${PCA}" > "${PQ}"
#   printf -- 'done\n'
#
#   if [[ ! -s ${PQ} ]]
#   then
#       error 'No variant in the PCA variants file is mentioned in the query VCF.'
#   fi
#
#   PR="${WORKDIR}/PR"
#   printf -- "${TEMPLATE}" "${PR}"
#   comm -1 -2 "${REFERENCE}" "${PCA}" > "${PR}"
#   printf -- 'done\n'
#
#   if [[ ! -s ${PR} ]]
#   then
#       error 'No variant in the PCA variants file is mentioned in the reference VCF.'
#   fi
#
#   NPR=$( wc --lines < "${PR}" )
#   NPCS=20
#   if (( NPR < 2 * NPCS + 1 ))
#   then
#
#       # NB: The flashpca program, invoked by PCATasks.PerformPCA, will fail if
#       #
#       #     min(NVARIANTS, NSAMPLES) < 2 * NPCS + 1
#       #
#       # ...where NVARIANTS and NSAMPLES are the numbers of variants (rows) and
#       # samples (data columns) in the reference VCF, and NPCS is the value of
#       # the command's --ndim flag (currently hard-coded as 20; see
#       # PCATasks.PerformPCA).  (The inequality above is derived from line 623 of
#       # https://github.com/gabraham/flashpca/blob/b8044f13607a072125828547684fde8b081d6191/flashpca.cpp .)
#       # In particular, flashpca will fail when NPR = NVARIANTS is less than
#       # 2 * NPCS + 1.  Hence the error condition detected here.
#
#       error 'The reference VCF mentions too few PCA variants for the PCA calculation.'
#   fi
#
#   NIXQ="${WORKDIR}/NIXQ"
#   printf -- "${TEMPLATE}" "${NIXQ}"
#   comm -2 -3 "${PQ}" "${PR}" > "${NIXQ}"
#   printf -- 'done\n'
#
#   NIXR="${WORKDIR}/NIXR"
#   printf -- "${TEMPLATE}" "${NIXR}"
#   comm -2 -3 "${PR}" "${PQ}" > "${NIXR}"
#   printf -- 'done\n'
#
#   NIX="${WORKDIR}/NIX"
#   printf -- "${TEMPLATE}" "${NIX}"
#   sort "${NIXQ}" "${NIXR}" > "${NIX}"
#   printf -- 'done\n'
#
#   if [[ ! -s ${NIX} ]]
#   then
#       exit 0
#   fi
#
#   (
#       sortvariants() {
#           local variants="${1}"
#           sort                    \
#               --field-separator=: \
#               --key=1,1V          \
#               --key=2,2n          \
#               --key=3,3           \
#               --key=4,4           \
#               "${variants}"
#       }
#
#       SEPARATOR=''
#       maybewarn() {
#           local nixx="${1}"
#           local label0="${2}"
#           local label1
#           if [[ ${label0} == query ]]
#           then
#               label1=reference
#           else
#               label1=query
#           fi
#           if [[ -s ${nixx} ]]
#           then
#               printf -- "${SEPARATOR}"
#               if [[ -z "${SEPARATOR}" ]]
#               then
#                   SEPARATOR=$'\n'
#               fi
#               printf -- 'WARNING: '
#               printf -- 'the following PCA variants from the '
#               printf -- '%s VCF do not appear in the %s VCF, '  \
#                         "${label0}" "${label1}"
#               printf -- 'and will be removed from the '
#               printf -- 'analysis.\n'
#               sortvariants "${nixx}"
#           fi
#       }
#
#       exec >>'~{WARNINGS}'
#       maybewarn "${NIXQ}" query
#       maybewarn "${NIXR}" reference
#   )
#
#   # ---------------------------------------------------------------------------
#
#   TEMP_PCA="${WORKDIR}/TEMP_PCA"
#   printf -- "${TEMPLATE}" "${TEMP_PCA}"
#   comm -2 -3 "${PCA}" "${NIX}" > "${TEMP_PCA}"
#   printf -- 'done\n'
#
#   WANTED_PCA="${WORKDIR}/WANTED_PCA"
#   printf -- "${TEMPLATE}" "${WANTED_PCA}"
#   join                                                \
#       -t $'\t'                                        \
#       "${TEMP_PCA}"                                   \
#       <(
#          perl -lpe '$_ = qq($_\t$.)' "${UNSORTEDPCA}" \
#            | sort                                     \
#                  --field-separator=$'\t'              \
#                  --key=1,1
#        )                                              \
#     | sort                                            \
#           --field-separator=$'\t'                     \
#           --key=2,2n                                  \
#     | cut --fields=1                                  \
#     > "${WANTED_PCA}"
#   printf -- 'done\n'
#
#   mkdir --verbose --parents "$( dirname '~{KEPT_PCA_VARIANTS}' )"
#   cp --verbose "${WANTED_PCA}" '~{KEPT_PCA_VARIANTS}'
#
#   # ---------------------------------------------------------------------------
#
#   printf -- '\n\n## WORKDIR:\n'
#   find '~{WORKDIR}' -type f | xargs ls -ltr
#
#   (
#       exec >&2
#       printf -- '\n\n### LINE COUNTS:\n'
#       find '~{WORKDIR}' -type f | xargs ls -1tr | xargs wc --lines
#   )
#
#   # ---------------------------------------------------------------------------
#
#   if ~{if debug then "true" else "false"}
#   then
#       tar -cvzf '~{ARCHIVE}' '~{WORKDIR}'
#   fi
#
#   # ---------------------------------------------------------------------------
#
#   printf -- 'FINAL STORAGE UTILIZATION:\n'
#   df --human
#
#   # ---------------------------------------------------------------------------
#
#   if ~{if !debug then "true" else "false"}
#   then
#       rm -rf '~{WORKDIR}'
#   fi
#
#   # ---------------------------------------------------------------------------
#   >>>
#
#   output {
#     File?        kept_pca_variants    = KEPT_PCA_VARIANTS
#     File?        warnings             = WARNINGS
#
#     Array[File?] workfiles            = WORKFILES
#     File?        workfiles_tgz        = ARCHIVE
#   }
#
#   runtime {
#     docker: "ubuntu:21.10"
#     disks : "local-disk ~{storage} HDD"
#   }
# }

task PurgeTmp {
  input {
    String tmp
    String sequencing_dummy = "ignored"
    String workspace
  }

  command <<<
  export WORKSPACE='~{workspace}'
  TMP="$( mapurl.sh '~{tmp}' )"
  rclone purge "${TMP}"
  >>>
  runtime {
    docker: "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/sharding:0.0.1"
  }
}
