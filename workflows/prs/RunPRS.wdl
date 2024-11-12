version 1.0

import "ScoringPart.wdl"
import "ScoringTasks.wdl"
import "HelperTasks.wdl"
import "Structs.wdl"

workflow RunPRSWorkflow {

  input { 
    File             query_vcf

    File             weights
    File             pca_variants
    File             reference_vcf
    Boolean          make_model_only = false                     

    String           name
  }

  call HelperTasks.GetBaseMemory as GetMemoryForQuery {
    input:
        vcf = query_vcf
  }

  call HelperTasks.GetBaseMemory as GetMemoryForReference {
    input:
        vcf = reference_vcf
  }

  call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInWeights {
    input:
        tsv        = weights
      , skipheader = true
  }

  call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInPcaVariants {
    input:
        tsv        = pca_variants
      , skipheader = false
  }

  call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
    input:
      vcf = query_vcf
  }

  call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInReferenceVcf {
    input:
      vcf = reference_vcf
  }

  call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
    input:
        vcf = RenameChromosomesInQueryVcf.renamed
      , mem = GetMemoryForQuery.gigabytes
  }

  call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
    input:
        vcf = RenameChromosomesInReferenceVcf.renamed
      , mem = GetMemoryForReference.gigabytes
  }

  call GetRegions {
    input:
        weights      = RenameChromosomesInWeights.renamed
      , pca_variants = RenameChromosomesInPcaVariants.renamed
      , query        = ExtractQueryVariants.ids
      , reference    = ExtractReferenceVariants.ids
  }

  call HelperTasks.GetBaseMemory {
    input:
      nvariants = GetRegions.nvariants
  }

  if (defined(GetRegions.query_regions)) {
    call HelperTasks.SubsetVcf as SubsetQueryVcf {
      input:
          inputvcf = RenameChromosomesInQueryVcf.renamed
        , regions  = select_first([GetRegions.query_regions])
        , label    = "query"
    }
  }

  if (defined(GetRegions.reference_regions)) {
    call HelperTasks.SubsetVcf as SubsetReferenceVcf {
      input:
          inputvcf = RenameChromosomesInReferenceVcf.renamed
        , regions  = select_first([GetRegions.reference_regions])
        , label    = "reference"
    }
  }

  WeightSet weight_set = object {
    linear_weights : RenameChromosomesInWeights.renamed
  }

  NamedWeightSet named_weight_set = object {
      condition_name : name
    , weight_set     : weight_set
  }

  Int base_memory  = GetBaseMemory.gigabytes
  Int pca_memory   = base_memory
  Int plink_memory = base_memory

  call ScoringPart.ScoringImputedDataset as ScoreQueryVcf {
    input:

        named_weight_set                     = named_weight_set
      , pruning_sites_for_pca                = RenameChromosomesInPcaVariants.renamed
      , imputed_array_vcf                    = select_first([SubsetQueryVcf.result,
                                                             RenameChromosomesInQueryVcf.renamed])
      , population_vcf                       = select_first([SubsetReferenceVcf.result,
                                                             RenameChromosomesInReferenceVcf.renamed])
      , make_model_only                      = make_model_only
      , redoPCA                              = true
      , adjustScores                         = true

      , basename                             = name
      , population_basename                  = "1kg"

      , extract_ids_plink_mem                = plink_memory
      , extract_ids_population_mem           = plink_memory
      , vcf_to_plink_mem                     = plink_memory
      , population_vcf_to_plink_mem          = plink_memory
      , scoring_mem                          = plink_memory
      , population_scoring_mem               = plink_memory
      , pca_memory                           = pca_memory
      , project_array_memory                 = pca_memory

      , population_loadings                  = "PLACEHOLDER__REQUIRED_BUT_NOT_USED"
      , population_pcs                       = "PLACEHOLDER__REQUIRED_BUT_NOT_USED"
      , population_meansd                    = "PLACEHOLDER__REQUIRED_BUT_NOT_USED"
  }

  output {
    File?   raw_scores                 = ScoreQueryVcf.raw_scores
    File?   adjusted_array_scores      = ScoreQueryVcf.adjusted_array_scores
    File    pc_projection              = select_first([ScoreQueryVcf.pc_projection])
    File    pc_plot                    = select_first([ScoreQueryVcf.pc_plot])

    # -------------------------------------------------------------------------

    File    model_parameters           = select_first([ScoreQueryVcf.model_parameters])
    File    training_variants          = select_first([ScoreQueryVcf.training_variants])
    # .........................................................................
    File    pcs                        = select_first([ScoreQueryVcf.pcs])
    File    pcloadings                 = select_first([ScoreQueryVcf.pcloadings])
    File    pcmeansd                   = select_first([ScoreQueryVcf.pcmeansd])

    # -------------------------------------------------------------------------

    File?   query_regions              = GetRegions.query_regions
    Boolean converged                  = select_first([ScoreQueryVcf.fit_converged])
    File    adjusted_population_scores = select_first([ScoreQueryVcf.adjusted_population_scores])

    # Int?    n_missing_sites_from_training = CompareScoredSitesToSitesUsedInTraining.n_missing_sites
    # File?   missing_sites_shifted_scores  = CombineMissingSitesAdjustedScores.missing_sites_shifted_scores
  }
}

# -------------------------------------------------------------------------------

task GetRegions {
  input {
    File    weights
    File    pca_variants
    File    query
    File    reference
    Boolean debug        = false
  }

  String        WORKDIR           = "WORK"
  String        ARCHIVE           = WORKDIR + ".tgz"
  String        OUTPUTDIR         = "OUTPUT"
  String        QUERY_REGIONS     = OUTPUTDIR + "/query_regions.tsv"
  String        REFERENCE_REGIONS = OUTPUTDIR + "/reference_regions.tsv"
  String        NVARIANTS         = OUTPUTDIR + "/nvariants"
  String        WARNINGS          = OUTPUTDIR + "/WARNINGS"
  Array[String] WORKFILES         = [
                                        WORKDIR + "/WEIGHTS"
                                      , WORKDIR + "/PCA"
                                      , WORKDIR + "/QUERY"
                                      , WORKDIR + "/REFERENCE"
                                      , WORKDIR + "/PQ"
                                      , WORKDIR + "/PR"
                                      , WORKDIR + "/NIXQ"
                                      , WORKDIR + "/NIXR"
                                      , WORKDIR + "/NIX"
                                      , WORKDIR + "/WQ"
                                      , WORKDIR + "/WR"
                                      , WORKDIR + "/WQ_U_WR"
                                      , WORKDIR + "/EXTRA"
                                      , WORKDIR + "/PQR"
                                      , WORKDIR + "/WANTED"
                                      , WORKDIR + "/QS"
                                      , WORKDIR + "/RS"
                                      , WORKDIR + "/QSP"
                                      , WORKDIR + "/RSP"
                                      , WORKDIR + "/REGIONS"
                                    ]
  Int           multiplier        = if debug then 20 else 10
  Int           storage           = 20 + multiplier * ceil(  size(weights     , "GB")
                                                           + size(pca_variants, "GB")
                                                           + size(query       , "GB")
                                                           + size(reference   , "GB"))

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  error() {
      local message="${1}"
      printf -- '%s\n' "${message}" >&2
      exit 1
  }

  # ---------------------------------------------------------------------------

  printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
  printf -- 'INITIAL STORAGE UTILIZATION:\n'
  df --human
  printf -- '\n'

  # ---------------------------------------------------------------------------

  WORKDIR='~{WORKDIR}'
  OUTPUTDIR='~{OUTPUTDIR}'
  UNSORTEDWEIGHTS='~{weights}'     # weights file (with headers row)
  UNSORTEDPCA='~{pca_variants}'    # pca_variants file (no header row)
  UNSORTEDREFERENCE='~{reference}' # CHROM:POS:REF:ALT from reference vcf
  UNSORTEDQUERY='~{query}'         # CHROM:POS:REF:ALT from query vcf

  # ---------------------------------------------------------------------------

  mkdir --verbose --parents "${WORKDIR}"
  mkdir --verbose --parents "${OUTPUTDIR}"

  WIDTH=$(( ${#WORKDIR} + 10 ))
  TEMPLATE="Creating %-${WIDTH}s ... "
  unset WIDTH

  WEIGHTS="${WORKDIR}/WEIGHTS"
  PCA="${WORKDIR}/PCA"
  REFERENCE="${WORKDIR}/REFERENCE"
  QUERY="${WORKDIR}/QUERY"

  printf -- "${TEMPLATE}" "${WEIGHTS}"
  # The perl segment of the following pipeline prints the first line only if
  # it ends in a numeric expression (i.e. the file does not have a headers row)
  cut --fields=1 "${UNSORTEDWEIGHTS}"       \
    | perl -lne '
        BEGIN {
          $::WEIGHT = qr/
             \b-?
             (?:\d+|\d*\.\d+|\d+\.\d*)
             (?:[eE][+-]?\d+)?\b
          /x;
        }
        print if $. > 1 || /$::WEIGHT\s*$/' \
    | sort --unique                         \
    > "${WEIGHTS}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${PCA}"
  cut --fields=1 "${UNSORTEDPCA}" | sort --unique > "${PCA}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${REFERENCE}"
  sort --unique "${UNSORTEDREFERENCE}" > "${REFERENCE}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${QUERY}"
  sort --unique "${UNSORTEDQUERY}" > "${QUERY}"
  printf -- 'done\n'

  # ---------------------------------------------------------------------------

  PQ="${WORKDIR}/PQ"
  printf -- "${TEMPLATE}" "${PQ}"
  comm -1 -2 "${QUERY}" "${PCA}" > "${PQ}"
  printf -- 'done\n'

  if [[ ! -s ${PQ} ]]
  then
      error 'No variant in the PCA variants file is mentioned in the query VCF.'
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
  comm -2 -3 "${PQ}" "${PR}" > "${NIXQ}"
  printf -- 'done\n'

  NIXR="${WORKDIR}/NIXR"
  printf -- "${TEMPLATE}" "${NIXR}"
  comm -2 -3 "${PR}" "${PQ}" > "${NIXR}"
  printf -- 'done\n'

  # if ! [[ -s ${NIXQ} ]] && ! [[ -s ${NIXR} ]]
  # then
  #     exit 0
  # fi

  NIX="${WORKDIR}/NIX"
  printf -- "${TEMPLATE}" "${NIX}"
  sort "${NIXQ}" "${NIXR}" > "${NIX}"
  printf -- 'done\n'

  if [[ -s ${NIX} ]]
  then

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
  fi

  # ---------------------------------------------------------------------------

  WR="${WORKDIR}/WR"
  WQ="${WORKDIR}/WQ"
  WQ_U_WR="${WORKDIR}/WQ_U_WR"
  EXTRA="${WORKDIR}/EXTRA"
  PQR="${WORKDIR}/PQR"
  WANTED="${WORKDIR}/WANTED"

  printf -- "${TEMPLATE}" "${WR}"
  comm -1 -2 "${REFERENCE}" "${WEIGHTS}" > "${WR}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${WQ}"
  comm -1 -2 "${QUERY}" "${WEIGHTS}" > "${WQ}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${WQ_U_WR}"
  sort --unique "${WQ}" "${WR}" > "${WQ_U_WR}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${EXTRA}"
  comm -2 -3 "${WQ_U_WR}" "${NIX}" > "${EXTRA}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${PQR}"
  comm -1 -2 "${PQ}" "${PR}" > "${PQR}"
  printf -- 'done\n'

  if [[ ! -s ${PQR} ]]
  then
      error 'No variant in the PCA variants file is mentioned in both the query and reference VCFs'
  fi

  printf -- "${TEMPLATE}" "${WANTED}"
  sort --unique "${PQR}" "${EXTRA}" > "${WANTED}"
  printf -- 'done\n'

  # ---------------------------------------------------------------------------

  QS="${WORKDIR}/QS"
  RS="${WORKDIR}/RS"
  QSP="${WORKDIR}/QSP"
  RSP="${WORKDIR}/RSP"

  printf -- "${TEMPLATE}" "${QS}"
  comm -1 -2 "${QUERY}" "${WANTED}" > "${QS}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${RS}"
  comm -1 -2 "${REFERENCE}" "${WANTED}" > "${RS}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${QSP}"
  comm -1 -2 "${PCA}" "${QS}" > "${QSP}"
  printf -- 'done\n'
  printf -- "${TEMPLATE}" "${RSP}"
  comm -1 -2 "${PCA}" "${RS}" > "${RSP}"
  printf -- 'done\n'

  if ! diff "${QSP}" "${RSP}" > /dev/null
  then
      (
          exec >&2
          printf -- 'INTERNAL ERROR: UNEXPECTED MISMATCHES:'
          diff "${QSP}" "${RSP}" || true
      )
      exit 1
  fi

  # ----------------------------------------------------------------------------

  REGIONS="${WORKDIR}/REGIONS"

  printf -- "${TEMPLATE}" "${REGIONS}"

  tr : $'\t' < "${WANTED}"        \
    | cut --fields=1,2            \
    | sort                        \
          --field-separator=$'\t' \
          --key=1,1V              \
          --key=2,2n              \
    > "${REGIONS}"
  printf -- 'done\n'

  # if [[ -s ${NIXQ} ]]
  # then
  #     cp "${REGIONS}" '~{QUERY_REGIONS}'
  # fi
  cp "${REGIONS}" '~{QUERY_REGIONS}'

  if [[ -s ${NIXR} ]]
  then
      cp "${REGIONS}" '~{REFERENCE_REGIONS}'
  fi

  printf -- '\n\n## WORKDIR:\n'
  find '~{WORKDIR}' -type f | xargs ls -ltr

  if ~{if debug then "true" else "false"}
  then
      tar -cvzf '~{ARCHIVE}' '~{WORKDIR}'
      (
          exec >&2
          printf -- '\n\n### LINE COUNTS:\n'
          find '~{WORKDIR}' -type f | xargs ls -1tr | xargs wc --lines
      )
  fi

  # ---------------------------------------------------------------------------

  wc --lines < "${REGIONS}" > '~{NVARIANTS}'

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
    # File? QUERY_REGIONS
    File         query_regions     = QUERY_REGIONS
    File?        reference_regions = REFERENCE_REGIONS
    File?        warnings          = WARNINGS
    Int          nvariants         = read_int(NVARIANTS)

    Array[File?] workfiles         = WORKFILES
    File?        workfiles_tgz     = ARCHIVE
  }

  runtime {
    docker: "ubuntu:21.10"
    disks : "local-disk ~{storage} HDD"
  }
}
