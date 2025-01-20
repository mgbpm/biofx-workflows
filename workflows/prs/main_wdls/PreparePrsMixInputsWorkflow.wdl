version 1.0

import "../tasks/HelperTasks.wdl"
import "../palantir/ScoringTasks.wdl"
import "../subwdls/BuildReferenceVcfWorkflow.wdl"

workflow PreparePrsMixInputs {
  input {
    Array[File] weights
    File        pca_variants
    String      workspace
    String      shards_root
    String      work_bucket
    Boolean     nocleanup        = false
    Int         nbatches         = 500
    Array[File] query_vcfs
    String      prs_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250122"
  }

  call BuildReferenceVcfWorkflow.BuildReferenceVcf {
    input:
        weights      = weights
      , pca_variants = pca_variants
      , workspace    = workspace
      , shards_root  = shards_root
      , work_bucket  = work_bucket
      , nocleanup    = nocleanup
      , nbatches     = nbatches
  }

  call HelperTasks.GetBaseMemory as GetMemoryForReference {
    input:
        vcf = BuildReferenceVcf.reference_vcf
  }

  call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
    input:
        vcf = BuildReferenceVcf.reference_vcf
      , mem = GetMemoryForReference.gigabytes
  }

  scatter (query_vcf in query_vcfs) {
    call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
      input:
          vcf = query_vcf
    }

    call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
      input:
          vcf = query_vcf
    }

    call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
      input:
          vcf = RenameChromosomesInQueryVcf.renamed
        , mem = GetMemoryForQueryFromVcf.gigabytes
    }
  }

  call HelperTasks.GetTotalSize as FootprintOfVariantFiles {
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
        pca_variants = BuildReferenceVcf.renamed_pca_variants
      , reference    = ExtractReferenceVariants.ids
      , union        = Union.result
      , intersection = Intersection.result
  }

  output {
    File        kept_pca_variants  = select_first([MaybeTrimPcaVariants.kept_pca_variants,
                                                   BuildReferenceVcf.renamed_pca_variants])
    Array[File] renamed_weights    = BuildReferenceVcf.renamed_weights
    Array[File] renamed_query_vcfs = RenameChromosomesInQueryVcf.renamed
    File        reference_vcf      = BuildReferenceVcf.reference_vcf
    File        reference_tbi      = BuildReferenceVcf.reference_tbi
    File        regions            = BuildReferenceVcf.regions
  }
}

# -----------------------------------------------------------------------------

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
