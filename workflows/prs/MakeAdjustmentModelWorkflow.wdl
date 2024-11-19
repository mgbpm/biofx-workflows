version 1.0

import "PCATasks.wdl"
import "ScoringTasks.wdl"
import "TrainAncestryAdjustmentModel.wdl"
import "Structs.wdl"
import "HelperTasks.wdl"

workflow MakeAdjustmentModel {
  input {
    File   weights
    File   pca_variants
    File   reference_vcf
    File   query_file
    String name
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

  call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInReferenceVcf {
    input:
        vcf = reference_vcf
  }

  call HelperTasks.GetBaseMemory as GetMemoryForReference {
    input:
        vcf = reference_vcf
  }

  call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
    input:
        vcf = RenameChromosomesInReferenceVcf.renamed
      , mem = GetMemoryForReference.gigabytes
  }

  Boolean isvcf = basename(query_file) != basename(query_file, ".vcf.gz")

  if (isvcf) {

      call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
        input:
            vcf = query_file
      }

      call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
        input:
            vcf = query_file
      }

      call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
        input:
            vcf = RenameChromosomesInQueryVcf.renamed
          , mem = GetMemoryForQueryFromVcf.gigabytes
      }
  }

  if (!isvcf) {

      call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInQueryVariants {
        input:
            tsv        = query_file
          , skipheader = false
      }
  }

  call GetRegions {
    input:
        weights      = RenameChromosomesInWeights.renamed
      , pca_variants = RenameChromosomesInPcaVariants.renamed
      , query        = select_first([ExtractQueryVariants.ids,
                                     RenameChromosomesInQueryVariants.renamed])
      , reference    = ExtractReferenceVariants.ids
  }

  call HelperTasks.GetBaseMemory {
    input:
        nvariants = GetRegions.n_reference_variants
  }

  call HelperTasks.SubsetVcf as SubsetReferenceVcf {
    input:
        inputvcf = RenameChromosomesInReferenceVcf.renamed
      , regions  = GetRegions.reference_regions
      , label    = "reference"
  }

  String reference_basename = basename(reference_vcf, ".vcf.gz")

  call PCATasks.ArrayVcfToPlinkDataset as ReferenceBed {
    input:
        vcf             = select_first([SubsetReferenceVcf.result,
                                        RenameChromosomesInReferenceVcf.renamed])
      , pruning_sites   = RenameChromosomesInPcaVariants.renamed
      , basename        = reference_basename
      , mem             = GetBaseMemory.gigabytes
      , subset_to_sites = GetRegions.query_variants
  }

  call PCATasks.PerformPCA {
    input:
        bed      = ReferenceBed.bed
      , bim      = ReferenceBed.bim
      , fam      = ReferenceBed.fam
      , basename = reference_basename
      , mem      = GetMemoryForReference.gigabytes
  }

  WeightSet weight_set = object {
    linear_weights : RenameChromosomesInWeights.renamed
  }

  NamedWeightSet named_weight_set = object {
    condition_name : name,
    weight_set     : weight_set
  }

  call TrainAncestryAdjustmentModel.TrainAncestryAdjustmentModel as TrainModel {
    input:
        named_weight_set    = named_weight_set
      , population_pcs      = PerformPCA.pcs
      , population_vcf      = RenameChromosomesInReferenceVcf.renamed
      , population_basename = basename(reference_vcf, ".vcf.gz")
      , sites               = GetRegions.query_variants
  }

  call BundleAdjustmentModel {
    input:
        model_data = object {
            parameters           : TrainModel.fitted_params
          , training_variants    : TrainModel.sites_used_in_scoring

          , principal_components : PerformPCA.pcs
          , loadings             : PerformPCA.pc_loadings
          , meansd               : PerformPCA.mean_sd

          , weights              : RenameChromosomesInWeights.renamed
          , pca_variants         : RenameChromosomesInPcaVariants.renamed

          , base_memory          : GetMemoryForReference.gigabytes
          , query_regions        : GetRegions.query_regions
        }
  }

  output {
    File    adjustment_model_manifest = BundleAdjustmentModel.manifest
    Boolean converged                 = TrainModel.fit_converged
    File    raw_reference_scores      = TrainModel.raw_population_scores
    File    adjusted_reference_scores = TrainModel.adjusted_population_scores
  }
}

# -----------------------------------------------------------------------------

task GetRegions {
  input {
    File    weights
    File    pca_variants
    File    query
    File    reference
    Boolean debug        = false
  }

  String        WORKDIR                 = "WORK"
  String        ARCHIVE                 = WORKDIR + ".tgz"
  String        OUTPUTDIR               = "OUTPUT"
  String        QUERY_REGIONS           = OUTPUTDIR + "/query_regions.tsv"
  String        REFERENCE_REGIONS       = OUTPUTDIR + "/reference_regions.tsv"
  String        KEPT_QUERY_VARIANTS     = OUTPUTDIR + "/kept_query_variants.txt"
  String        KEPT_REFERENCE_VARIANTS = OUTPUTDIR + "/kept_query_variants.txt"
  String        N_QUERY_VARIANTS        = OUTPUTDIR + "/n_query_variants"
  String        N_REFERENCE_VARIANTS    = OUTPUTDIR + "/n_reference_variants"
  String        WARNINGS                = OUTPUTDIR + "/WARNINGS"
  Array[String] WORKFILES               = [
                                              WORKDIR + "/WEIGHTS"
                                            , WORKDIR + "/PCA"
                                            , WORKDIR + "/QUERY"
                                            , WORKDIR + "/REFERENCE"
                                            , WORKDIR + "/PQ"
                                            , WORKDIR + "/PR"
                                            , WORKDIR + "/NIXQ"
                                            , WORKDIR + "/NIXR"
                                            , WORKDIR + "/WQ"
                                            , WORKDIR + "/WR"
                                            , WORKDIR + "/EXTRA_FOR_Q"
                                            , WORKDIR + "/EXTRA_FOR_R"
                                            , WORKDIR + "/PQR"
                                            , WORKDIR + "/WANTED_Q"
                                            , WORKDIR + "/WANTED_R"
                                            , WORKDIR + "/QS"
                                            , WORKDIR + "/RS"
                                            , WORKDIR + "/QSP"
                                            , WORKDIR + "/RSP"
                                          ]
  Int           multiplier          = if debug then 20 else 10
  Int           storage             = 30 + multiplier * ceil(  size(weights     , "GB")
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

  if [[ -s ${NIXQ} ]] || [[ -s ${NIXR} ]]
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
  EXTRA_FOR_Q="${WORKDIR}/EXTRA_FOR_Q"
  EXTRA_FOR_R="${WORKDIR}/EXTRA_FOR_R"
  PQR="${WORKDIR}/PQR"
  WANTED_Q="${WORKDIR}/WANTED_Q"
  WANTED_R="${WORKDIR}/WANTED_R"

  printf -- "${TEMPLATE}" "${WR}"
  comm -1 -2 "${REFERENCE}" "${WEIGHTS}" > "${WR}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${WQ}"
  comm -1 -2 "${QUERY}" "${WEIGHTS}" > "${WQ}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${EXTRA_FOR_Q}"
  comm -2 -3 "${WQ}" "${NIXQ}" > "${EXTRA_FOR_Q}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${EXTRA_FOR_R}"
  comm -2 -3 "${WR}" "${NIXR}" > "${EXTRA_FOR_R}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${PQR}"
  comm -1 -2 "${PQ}" "${PR}" > "${PQR}"
  printf -- 'done\n'

  if [[ ! -s ${PQR} ]]
  then
      error 'No variant in the PCA variants file is mentioned in both the query and reference VCFs'
  fi

  printf -- "${TEMPLATE}" "${WANTED_Q}"
  sort --unique "${PQR}" "${EXTRA_FOR_Q}" > "${WANTED_Q}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${WANTED_R}"
  sort --unique "${PQR}" "${EXTRA_FOR_R}" > "${WANTED_R}"
  printf -- 'done\n'

  mkdir --verbose --parents "$( dirname '~{KEPT_QUERY_VARIANTS}' )"
  cp --verbose "${WANTED_Q}" '~{KEPT_QUERY_VARIANTS}'

  # ---------------------------------------------------------------------------

  QS="${WORKDIR}/QS"
  RS="${WORKDIR}/RS"
  QSP="${WORKDIR}/QSP"
  RSP="${WORKDIR}/RSP"

  printf -- "${TEMPLATE}" "${QS}"
  comm -1 -2 "${QUERY}" "${WANTED_Q}" > "${QS}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${RS}"
  comm -1 -2 "${REFERENCE}" "${WANTED_R}" > "${RS}"
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

  REGIONS="${WORKDIR}/REGIONS_Q"
  REGIONS="${WORKDIR}/REGIONS_R"

  printf -- "${TEMPLATE}" '~{QUERY_REGIONS}'
  tr : $'\t' < "${WANTED_Q}"      \
    | cut --fields=1,2            \
    | sort                        \
          --field-separator=$'\t' \
          --key=1,1V              \
          --key=2,2n              \
    > '~{QUERY_REGIONS}'
  printf -- 'done\n'

  printf -- "${TEMPLATE}" '~{REFERENCE_REGIONS}'
  tr : $'\t' < "${WANTED_R}"      \
    | cut --fields=1,2            \
    | sort                        \
          --field-separator=$'\t' \
          --key=1,1V              \
          --key=2,2n              \
    > '~{REFERENCE_REGIONS}'
  printf -- 'done\n'

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

  wc --lines < "${WANTED_Q}" > '~{N_QUERY_VARIANTS}'
  wc --lines < "${WANTED_R}" > '~{N_REFERENCE_VARIANTS}'

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
    File         query_regions        = QUERY_REGIONS
    File         reference_regions    = REFERENCE_REGIONS
    File         query_variants       = KEPT_QUERY_VARIANTS
    File         reference_variants   = KEPT_REFERENCE_VARIANTS
    Int          n_query_variants     = read_int(N_QUERY_VARIANTS)
    Int          n_reference_variants = read_int(N_REFERENCE_VARIANTS)
    File?        warnings             = WARNINGS

    Array[File?] workfiles            = WORKFILES
    File?        workfiles_tgz        = ARCHIVE
  }

  runtime {
    docker: "ubuntu:21.10"
    disks : "local-disk ~{storage} HDD"
  }
}

task BundleAdjustmentModel {
  input {
    AdjustmentModelData model_data
  }

  command {}

  output {
    File manifest = write_json(model_data)
  }

  runtime {
    docker: "ubuntu:21.10"
  }
}

# -----------------------------------------------------------------------------
