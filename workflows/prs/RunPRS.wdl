version 1.0

import "ScoringPart.wdl"
import "Structs.wdl"

workflow RunPRSWorkflow {

  input { 
    File             weights
    File             pca_variants
    File             imputed_vcf
    File             reference_vcf
                     
    String           name

    Map[String, Int] memory_specs = {}
  }

  call Baseline

  call ExtractVariants as ExtractImputedVariants {
    input:
      vcf = imputed_vcf
  }

  call ExtractVariants as ExtractReferenceVariants {
    input:
      vcf = reference_vcf
  }

  call GetRegions {
    input:
      weights      = weights,
      pca_variants = pca_variants,
      imputed      = ExtractImputedVariants.variants,
      reference    = ExtractReferenceVariants.variants
  }

  if (defined(GetRegions.imputed_regions)) {
    call SubsetVcf as SubsetImputedVcf {
      input:
        inputvcf = imputed_vcf,
        regions  = select_first([GetRegions.imputed_regions]),
        label    = "imputed"
    }
  }

  if (defined(GetRegions.reference_regions)) {
    call SubsetVcf as SubsetReferenceVcf {
      input:
        inputvcf = reference_vcf,
        regions  = select_first([GetRegions.reference_regions]),
        label    = "reference"
    }
  }

  WeightSet weight_set = object {
    linear_weights : weights
  }

  NamedWeightSet named_weight_set = object {
    condition_name : name,
    weight_set     : weight_set
  }

  # Map[String, Int] memory_specs_default = {
  #                                           "i_extract"      : 32,
  #                                           "i_scoring"      : 16,
  #                                           "i_vcf_to_plink" : 16,
  #                                           "r_extract"      : 64,
  #                                           "r_scoring"      : 64,
  #                                           "r_vcf_to_plink" : 64
  #                                         }

  call ScoringPart.ScoringImputedDataset {
    input:

      named_weight_set                     = named_weight_set,
      pruning_sites_for_pca                = pca_variants,
      imputed_array_vcf                    = select_first([SubsetImputedVcf.result  , imputed_vcf  ]),
      population_vcf                       = select_first([SubsetReferenceVcf.result, reference_vcf]),
      redoPCA                              = true,
      adjustScores                         = true,

      basename                             = name,
      population_basename                  = "1kg",

      # ExtractIDsPlink.mem                  = select_first([memory_specs        ["i_extract"],
      #                                                      memory_specs_default["i_extract"]])
      # scoring_mem                          = 16,
      # vcf_to_plink_mem                     = 16,
      # 
      # ExtractIDsPopulation.mem             = 64,
      # population_scoring_mem               = 64,
      # PopulationArrayVcfToPlinkDataset.mem = 64

      population_loadings                  = "PLACEHOLDER__REQUIRED_BY_BUGGY_CODE",
      population_pcs                       = "PLACEHOLDER__REQUIRED_BY_BUGGY_CODE",
      population_meansd                    = "PLACEHOLDER__REQUIRED_BY_BUGGY_CODE"
  }

  output {
    File    raw_scores                    = ScoringImputedDataset.raw_scores
    Boolean fit_converged                 = select_first([ScoringImputedDataset.fit_converged])
    File    adjusted_array_scores         = select_first([ScoringImputedDataset.adjusted_array_scores])
    File    adjusted_population_scores    = select_first([ScoringImputedDataset.adjusted_population_scores])
    File    pc_projection                 = select_first([ScoringImputedDataset.pc_projection])
    File    pc_plot                       = select_first([ScoringImputedDataset.pc_plot])

    # Int?    n_missing_sites_from_training = CompareScoredSitesToSitesUsedInTraining.n_missing_sites
    # File?   missing_sites_shifted_scores  = CombineMissingSitesAdjustedScores.missing_sites_shifted_scores
  }
}

# -------------------------------------------------------------------------------

task Baseline {

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  printf -- 'ENVIRONMENT:\n'
  typeset -p
  printf -- '\n'

  printf -- 'INITIAL STORAGE UTILIZATION:\n'
  df --human

  # ---------------------------------------------------------------------------
  >>>

  runtime {
    docker: "ubuntu:21.10"
  }
}

# -------------------------------------------------------------------------------

task ExtractVariants {
  input {
    File vcf
  }

  String OUTPUTDIR = "OUTPUT"
  String OUTPUT    = OUTPUTDIR + "/variants.txt"
  Int    storage   = 20 + 2 * ceil(size(vcf, "GB"))

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
  printf -- 'INITIAL STORAGE UTILIZATION:\n'
  df --human
  printf -- '\n'

  # ---------------------------------------------------------------------------

  mkdir --verbose --parents '~{OUTPUTDIR}'

  zcat --force '~{vcf}'    \
    | grep                 \
          --invert-match   \
          --perl-regexp    \
          '^#'             \
    | cut --fields=1,2,4,5 \
    | tr $'\t' :           \
    > '~{OUTPUT}'

  # ---------------------------------------------------------------------------

  printf -- 'FINAL STORAGE UTILIZATION:\n'
  df --human

  # ---------------------------------------------------------------------------
  >>>

  output {
    File variants = OUTPUT
  }

  runtime {
    docker: "ubuntu:21.10"
    disks : "local-disk ~{storage} HDD"
  }
}


task GetRegions {
  input {
    File weights
    File pca_variants
    File imputed
    File reference
  }

  String OUTPUTDIR         = "OUTPUT"
  String IMPUTED_REGIONS   = OUTPUTDIR + "/imputed_regions.tsv"
  String REFERENCE_REGIONS = OUTPUTDIR + "/reference_regions.tsv"
  String WARNINGS          = OUTPUTDIR + "/WARNINGS"
  Int    storage           = 20 + 10 * ceil(  size(weights     , "GB")
                                            + size(pca_variants, "GB")
                                            + size(imputed     , "GB")
                                            + size(reference   , "GB"))

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
  printf -- 'INITIAL STORAGE UTILIZATION:\n'
  df --human
  printf -- '\n'

  # ---------------------------------------------------------------------------

  OUTPUTDIR='~{OUTPUTDIR}'
  UNSORTEDWEIGHTS='~{weights}'     # weights file (with headers row)
  UNSORTEDPCA='~{pca_variants}'    # pca_variants file (no header row)
  UNSORTEDREFERENCE='~{reference}' # CHROM:POS:REF:ALT from reference vcf
  UNSORTEDIMPUTED='~{imputed}'     # CHROM:POS:REF:ALT from imputed vcf

  # ---------------------------------------------------------------------------

  mkdir --verbose --parents "${OUTPUTDIR}"

  WORKDIR="$( mktemp --directory )"
  WIDTH=$(( ${#WORKDIR} + 10 ))
  TEMPLATE="Creating %-${WIDTH}s ... "
  unset WIDTH

  WEIGHTS="${WORKDIR}/WEIGHTS"
  PCA="${WORKDIR}/PCA"
  REFERENCE="${WORKDIR}/REFERENCE"
  IMPUTED="${WORKDIR}/IMPUTED"

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
  sort --unique  "${UNSORTEDREFERENCE}" > "${REFERENCE}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${IMPUTED}"
  sort --unique  "${UNSORTEDIMPUTED}"   > "${IMPUTED}"
  printf -- 'done\n'

  # ---------------------------------------------------------------------------

  PI="${WORKDIR}/PI"
  printf -- "${TEMPLATE}" "${PI}"
  comm -1 -2 "${IMPUTED}" "${PCA}" > "${PI}"
  printf -- 'done\n'

  PR="${WORKDIR}/PR"
  printf -- "${TEMPLATE}" "${PR}"
  comm -1 -2 "${REFERENCE}" "${PCA}" > "${PR}"
  printf -- 'done\n'

  NIXI="${WORKDIR}/NIXI"
  printf -- "${TEMPLATE}" "${NIXI}"
  comm -2 -3 "${PI}" "${PR}" > "${NIXI}"
  printf -- 'done\n'

  NIXR="${WORKDIR}/NIXR"
  printf -- "${TEMPLATE}" "${NIXR}"
  comm -2 -3 "${PR}" "${PI}" > "${NIXR}"
  printf -- 'done\n'

  # if ! [[ -s ${NIXI} ]] && ! [[ -s ${NIXR} ]]
  # then
  #     exit 0
  # fi

  NIX="${WORKDIR}/NIX"
  printf -- "${TEMPLATE}" "${NIX}"
  sort "${NIXI}" "${NIXR}" > "${NIX}"
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
              if [[ ${label0} == imputed ]]
              then
                  label1=reference
              else
                  label1=imputed
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
          maybewarn "${NIXI}" imputed
          maybewarn "${NIXR}" reference
      )
  fi

  # ---------------------------------------------------------------------------

  WR="${WORKDIR}/WR"
  WI="${WORKDIR}/WI"
  WR_U_WI="${WORKDIR}/WR_U_WI"
  EXTRA="${WORKDIR}/EXTRA"
  PRI="${WORKDIR}/PRI"
  WANTED="${WORKDIR}/WANTED"

  printf -- "${TEMPLATE}" "${WR}"
  comm -1 -2 "${REFERENCE}" "${WEIGHTS}" > "${WR}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${WI}"
  comm -1 -2 "${IMPUTED}"   "${WEIGHTS}" > "${WI}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${WR_U_WI}"
  sort --unique "${WR}" "${WI}" > "${WR_U_WI}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${EXTRA}"
  comm -2 -3 "${WR_U_WI}" "${NIX}" > "${EXTRA}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${PRI}"
  comm -1 -2 "${PR}" "${PI}" > "${PRI}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${WANTED}"
  sort --unique "${PRI}" "${EXTRA}" > "${WANTED}"
  printf -- 'done\n'

  # ---------------------------------------------------------------------------

  IS="${WORKDIR}/IS"
  RS="${WORKDIR}/RS"
  ISP="${WORKDIR}/ISP"
  RSP="${WORKDIR}/RSP"

  printf -- "${TEMPLATE}" "${IS}"
  comm -1 -2 "${IMPUTED}"   "${WANTED}" > "${IS}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${RS}"
  comm -1 -2 "${REFERENCE}" "${WANTED}" > "${RS}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${ISP}"
  comm -1 -2 "${PCA}" "${IS}" > "${ISP}"
  printf -- 'done\n'
  printf -- "${TEMPLATE}" "${RSP}"
  comm -1 -2 "${PCA}" "${RS}" > "${RSP}"
  printf -- 'done\n'

  if ! diff "${ISP}" "${RSP}" > /dev/null
  then
      (
          exec >&2
          printf -- 'INTERNAL ERROR: UNEXPECTED MISMATCHES:'
          diff "${ISP}" "${RSP}" || true
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

  # if [[ -s ${NIXI} ]]
  # then
  #     cp "${REGIONS}" '~{IMPUTED_REGIONS}'
  # fi
  cp "${REGIONS}" '~{IMPUTED_REGIONS}'

  if [[ -s ${NIXR} ]]
  then
      cp "${REGIONS}" '~{REFERENCE_REGIONS}'
  fi

  # ---------------------------------------------------------------------------

  printf -- 'FINAL STORAGE UTILIZATION:\n'
  df --human

  # ---------------------------------------------------------------------------
  >>>

  output {
    # File? IMPUTED_REGIONS
    File  imputed_regions   = IMPUTED_REGIONS
    File? reference_regions = REFERENCE_REGIONS
    File? warnings          = WARNINGS
  }

  runtime {
    docker: "ubuntu:21.10"
    disks : "local-disk ~{storage} HDD"
  }
}

task SubsetVcf {
  input {
    File    inputvcf
    File    regions
    String  label     = "data"
    Boolean nocleanup = false
  }

  String OUTPUTDIR = "OUTPUT"
  String OUTPUTVCF = OUTPUTDIR + "/" + label + ".vcf.gz"
  Int    storage   = 20 + 3 * ceil(size(inputvcf, "GB"))

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  set -o xtrace

  # ---------------------------------------------------------------------------

  printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
  printf -- 'INITIAL STORAGE UTILIZATION:\n'
  df --human
  printf -- '\n'

  # ---------------------------------------------------------------------------

  mkdir --verbose --parents '~{OUTPUTDIR}'
  WORKDIR="$( mktemp --directory )"
  INPUTVCF="${WORKDIR}/input.vcf.gz"

  ln --symbolic --verbose '~{inputvcf}' "${INPUTVCF}"

  bcftools index      \
      --force         \
      --tbi           \
      "${INPUTVCF}"

  cleanup() {

      bcftools                   \
          norm                   \
          --multiallelics -any   \
          --no-version           \
          --output-type v        \
    | bcftools                   \
          annotate               \
          --remove 'INFO,FORMAT' \
          --no-version           \
          --output-type v        \
    | perl -lne '
        BEGIN { $, = "\t"; }
        @::F = split( $,, $_, -1 );
        if (/^#/) {
          if (/^##contig=<ID=/) {
            unless ( $done ) {
              $done = 1;
              for $chromosome ( 1 .. 22, q(X) ) {
                print qq(##contig=<ID=$chromosome>);
              }
            }
          }
          else {
            print;
          }
        }
        else {
          $::F[ 2 ] = join(
            q(:),
            @::F[ 0, 1, 3, 4 ]
          );
          print @::F;
        }
      '
  }

  if ~{if nocleanup then "true" else "false"}
  then
      POSTPROCESS=cat
  else
      POSTPROCESS=cleanup
  fi

  bcftools                             \
      view                             \
      --no-version                     \
      --output-type v                  \
      --regions-file '~{regions}'      \
      "${INPUTVCF}"                    \
    | "${POSTPROCESS}"                 \
    | bcftools                         \
          view                         \
          --no-version                 \
          --output-type z              \
          --output-file '~{OUTPUTVCF}'

  # ---------------------------------------------------------------------------

  printf -- 'FINAL STORAGE UTILIZATION:\n'
  df --human

  # ---------------------------------------------------------------------------
  >>>

  output {
    File result = OUTPUTVCF
  }

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    disks : "local-disk ~{storage} HDD"
  }
}
