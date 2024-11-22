version 1.0

import "PCATasks.wdl"
import "ScoringTasks.wdl"
import "TrainAncestryAdjustmentModel.wdl"
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

  File query_variants = select_first([ExtractQueryVariants.ids,
                                      RenameChromosomesInQueryVariants.renamed])

  call MaybeTrimPcaVariants {
    input:
        pca_variants = RenameChromosomesInPcaVariants.renamed
      , reference    = ExtractReferenceVariants.ids
      , query        = query_variants
  }

  File   kept_pca_variants  = select_first([MaybeTrimPcaVariants.kept_pca_variants, 
                                            RenameChromosomesInPcaVariants.renamed])
  String reference_basename = basename(reference_vcf, ".vcf.gz")

  call PCATasks.ArrayVcfToPlinkDataset as ReferenceBed {
    input:
        vcf             = RenameChromosomesInReferenceVcf.renamed
      , pruning_sites   = kept_pca_variants
      , basename        = reference_basename
      , mem             = GetMemoryForReference.gigabytes
      , subset_to_sites = query_variants
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
      , population_basename = reference_basename
      , sites               = query_variants
  }

  call BundleAdjustmentModel {

    # NB: All the values to which "" is being prepended in the object
    # below are of File type.  Prepending "" to these values
    # effectively converts them into values of String type.  This
    # conversion is necessary to prevent Cromwell from localizing the
    # corresponding files.  In this case, such localization would be a
    # wasteful operation, since the task does not need these files at
    # all; it needs only their locations.

    input:
        model_data = object {
            parameters            : "" + TrainModel.fitted_params
          , training_variants     : "" + TrainModel.sites_used_in_scoring

          , principal_components  : "" + PerformPCA.pcs
          , loadings              : "" + PerformPCA.pc_loadings
          , meansd                : "" + PerformPCA.mean_sd

          , weights               : "" + RenameChromosomesInWeights.renamed
          , pca_variants          : "" + kept_pca_variants
          , original_pca_variants : "" + pca_variants

          , base_memory           :      GetMemoryForReference.gigabytes
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

task MaybeTrimPcaVariants {
  input {
    File    pca_variants
    File    reference
    File    query
    Boolean debug        = false
  }

  String        WORKDIR                 = "WORK"
  String        ARCHIVE                 = WORKDIR + ".tgz"
  String        OUTPUTDIR               = "OUTPUT"
  String        KEPT_PCA_VARIANTS       = OUTPUTDIR + "/kept_pca_variants.txt"
  String        WARNINGS                = OUTPUTDIR + "/WARNINGS"
  Array[String] WORKFILES               = [
                                              WORKDIR + "/PCA"
                                            , WORKDIR + "/QUERY"
                                            , WORKDIR + "/REFERENCE"
                                            , WORKDIR + "/PQ"
                                            , WORKDIR + "/PR"
                                            , WORKDIR + "/NIXQ"
                                            , WORKDIR + "/NIXR"
                                            , WORKDIR + "/NIX"
                                            , WORKDIR + "/TEMP_PCA"
                                            , WORKDIR + "/WANTED_PCA"
                                          ]
  Int           multiplier          = if debug then 20 else 10
  Int           storage             = 30 + multiplier * ceil(  size(pca_variants, "GB")
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
  UNSORTEDPCA='~{pca_variants}'    # pca_variants file (no header row)
  UNSORTEDREFERENCE='~{reference}' # CHROM:POS:REF:ALT from reference vcf
  UNSORTEDQUERY='~{query}'         # CHROM:POS:REF:ALT from query vcf

  # ---------------------------------------------------------------------------

  mkdir --verbose --parents "${WORKDIR}"
  mkdir --verbose --parents "${OUTPUTDIR}"

  WIDTH=$(( ${#WORKDIR} + 10 ))
  TEMPLATE="Creating %-${WIDTH}s ... "
  unset WIDTH

  PCA="${WORKDIR}/PCA"
  REFERENCE="${WORKDIR}/REFERENCE"
  QUERY="${WORKDIR}/QUERY"

  printf -- "${TEMPLATE}" "${PCA}"
  sort --unique "${UNSORTEDPCA}" > "${PCA}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${QUERY}"
  sort --unique "${UNSORTEDQUERY}" > "${QUERY}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${REFERENCE}"
  sort --unique "${UNSORTEDREFERENCE}" > "${REFERENCE}"
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
  find '~{WORKDIR}' -type f | xargs ls -ltr

  (
      exec >&2
      printf -- '\n\n### LINE COUNTS:\n'
      find '~{WORKDIR}' -type f | xargs ls -1tr | xargs wc --lines
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
  }

  runtime {
    docker: "ubuntu:21.10"
    disks : "local-disk ~{storage} HDD"
  }
}

task BundleAdjustmentModel {
  input {
    Object model_data
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
