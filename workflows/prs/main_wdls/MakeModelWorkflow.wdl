version 1.0

import "../tasks/PCATasks.wdl" as PCATasks
import "../tasks/ScoringTasks.wdl" as ScoringTasks
import "../subwdls/MixScoreWorkflow.wdl"
import "../subwdls/TrainModelWorkflow.wdl"
import "../tasks/HelperTasks.wdl" as HelperTasks
import "../../../steps/Utilities.wdl"

workflow MakeModelWorkflow {
  input {
    String      condition_code
    Array[File] variant_weights
    File        pca_variants
    File        reference_vcf
    File        query_file
    File?       score_weights
    Boolean     norename            = false
    String      ubuntu_docker_image = "ubuntu:21.10"
  }

  if ((length(variant_weights) > 1) && !defined(score_weights)) {
    call Utilities.FailTask as MissingScoreWeights {
      input:
        error_message = "Missing score_weights file. It is required when multiple variant weights are specified."
    }
  }

  if ((length(variant_weights) == 1) || defined(score_weights)) {
    if (defined(score_weights)) {
      call HelperTasks.CheckInputWeightFiles {
        input:
          score_weights = select_first([score_weights]),
          variant_weights = variant_weights,
          docker_image = ubuntu_docker_image
      }
    }
  
    if (!norename) {
      scatter (i in range(length(variant_weights))) {
        call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInWeights {
          input:
            tsv = variant_weights[i],
            skipheader = true
        }
      }
  
      call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInPcaVariants {
        input:
          tsv = pca_variants,
          skipheader = false
      }
  
      call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInReferenceVcf {
        input:
          vcf = reference_vcf
      }
    }
  
    Array[File] variant_weights_ = select_first([RenameChromosomesInWeights.renamed, variant_weights])
    File pca_variants_ = select_first([RenameChromosomesInPcaVariants.renamed, pca_variants])
    File reference_vcf_ = select_first([RenameChromosomesInReferenceVcf.renamed, reference_vcf])
  
    call HelperTasks.GetBaseMemory as GetMemoryForReference {
      input:
        vcf = reference_vcf
    }
  
    call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
      input:
        vcf = reference_vcf_,
        mem = GetMemoryForReference.gigabytes
    }
  
    Boolean isvcf = basename(query_file) != basename(query_file, ".vcf.gz")
  
    if (isvcf) {
      call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
        input:
          vcf = query_file
      }
  
      if (!norename) {
        call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
          input:
            vcf = query_file
        }
      }
  
      File query_vcf = select_first([RenameChromosomesInQueryVcf.renamed, query_file])
  
      call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
        input:
          vcf = select_first([RenameChromosomesInQueryVcf.renamed, query_file]),
          mem = GetMemoryForQueryFromVcf.gigabytes
      }
    }
    if (! isvcf) {
      if (!norename) {
        call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInQueryVariants {
          input:
            tsv = query_file,
            skipheader = false
        }
      }
      File query_file_ = select_first([RenameChromosomesInQueryVariants.renamed, query_file])
    }
  
    File query_variants = select_first([ExtractQueryVariants.ids, query_file_])
  
    call TrimPCAVariants {
      input:
        pca_variants = pca_variants_,
        reference = ExtractReferenceVariants.ids,
        query = query_variants,
        docker_image = ubuntu_docker_image
    }
  
    File kept_pca_variants = select_first([TrimPCAVariants.kept_pca_variants, pca_variants_])
  
    String reference_basename = basename(reference_vcf, ".vcf.gz")
  
    call PCATasks.ArrayVcfToPlinkDataset as ReferenceBed {
      input:
        vcf = reference_vcf_,
        pruning_sites = kept_pca_variants,
        subset_to_sites = query_variants,
        basename = reference_basename,
        mem = GetMemoryForReference.gigabytes
    }
    call PCATasks.PerformPCA {
      input:
        bed = ReferenceBed.bed,
        bim = ReferenceBed.bim,
        fam = ReferenceBed.fam,
        basename = reference_basename,
        mem = GetMemoryForReference.gigabytes
    }
  
    # Train the model
    scatter (weights_file in variant_weights_) {
      call ScoringTasks.ScoreVcf as ScorePopulationVcf {
        input:
          vcf = reference_vcf_,
          basename = basename(weights_file, ".tsv"),
          weights = weights_file,
          base_mem = GetMemoryForReference.gigabytes,
          sites = query_variants,
          chromosome_encoding = "MT"
      }

      Object scoring_inputs_ = object{
        variant_weights : "" + weights_file,
        training_variants: "" + ScorePopulationVcf.sites_scored
      }
    }

    if (length(variant_weights_) > 1) {
      call MixScoreWorkflow.MixScoreWorkflow as GetMixScore {
        input:
          output_basename = condition_code + "_" + basename(reference_vcf_),
          raw_scores = ScorePopulationVcf.score,
          score_weights = select_first([score_weights])
      }
    }

    File population_scores = select_first([GetMixScore.prs_mix_raw_score, ScorePopulationVcf.score[0]])

    call ScoringTasks.TrainAncestryModel {
      input:
        population_pcs = PerformPCA.pcs,
        population_scores = population_scores,
        output_basename = condition_code
    }
  
    if (defined(score_weights)) {
      call BundleAdjustmentModel as BundleMixModel {
        input:
          model_data = object {
            condition_code : condition_code,
            parameters : "" + TrainAncestryModel.fitted_params,
            scoring_inputs : scoring_inputs_,
            principal_components : "" + PerformPCA.pcs,
            loadings : "" + PerformPCA.pc_loadings,
            meansd : "" + PerformPCA.mean_sd,
            variant_weights : "" + variant_weights_,
            score_weights : "" + score_weights,
            pca_variants : "" + kept_pca_variants,
            original_pca_variants : "" + pca_variants,
            query_file : "" + select_first([query_vcf, query_file_]),
            base_memory : GetMemoryForReference.gigabytes
          },
          docker_image = ubuntu_docker_image
      }
    }
  
    if (!defined(score_weights)) {
      call BundleAdjustmentModel as BundleModel {
        input:
          model_data = object {
            condition_code : condition_code,
            parameters : "" + TrainAncestryModel.fitted_params,
            scoring_inputs : scoring_inputs_,
            principal_components : "" + PerformPCA.pcs,
            loadings : "" + PerformPCA.pc_loadings,
            meansd : "" + PerformPCA.mean_sd,
            variant_weights : "" + variant_weights_,
            pca_variants : "" + kept_pca_variants,
            original_pca_variants : "" + pca_variants,
            query_file : "" + select_first([query_vcf, query_file_]),
            base_memory : GetMemoryForReference.gigabytes
          },
          docker_image = ubuntu_docker_image
      }
    }
  }

  output {
    File         adjustment_model_manifest = select_first([BundleMixModel.manifest, BundleModel.manifest])
    Array[File]? raw_reference_scores      = ScorePopulationVcf.score
    File?        mixed_reference_scores    = GetMixScore.prs_mix_raw_score
    File?        adjusted_reference_scores = TrainAncestryModel.adjusted_population_scores
  }
}

task TrimPCAVariants {
  input {
    File    pca_variants
    File    reference
    File    query
    Boolean debug = false
    String  docker_image
  }

  String        WORKDIR           = "WORK"
  String        ARCHIVE           = WORKDIR + ".tgz"
  String        OUTPUTDIR         = "OUTPUT"
  String        KEPT_PCA_VARIANTS = OUTPUTDIR + "/kept_pca_variants.txt"
  String        WARNINGS          = OUTPUTDIR + "/WARNINGS"
  Array[String] WORKFILES         = [
    WORKDIR + "/PCA",
    WORKDIR + "/QUERY",
    WORKDIR + "/REFERENCE",
    WORKDIR + "/PQ",
    WORKDIR + "/PR",
    WORKDIR + "/NIXQ",
    WORKDIR + "/NIXR",
    WORKDIR + "/NIX",
    WORKDIR + "/TEMP_PCA",
    WORKDIR + "/WANTED_PCA"
    ]
  Int           multiplier        = if debug then 20 else 10
  Int           file_size         = ceil(size(pca_variants, "GB") + size(query, "GB")+ size(reference, "GB"))
  Int           final_disk_size   = 30 + multiplier * file_size

  command <<<
    set -o pipefail
    set -o errexit
    set -o nounset
  
    error() {
      local message="${1}"
      printf -- '%s\n' "${message}" >&2
      exit 1
    }
  
    printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{final_disk_size}'
    printf -- 'INITIAL STORAGE UTILIZATION:\n'
    df --human
    printf -- '\n'
  
    WORKDIR='~{WORKDIR}'
    OUTPUTDIR='~{OUTPUTDIR}'
    UNSORTEDPCA='~{pca_variants}'  # pca_variants file (no header row)
    UNSORTEDREFERENCE='~{reference}' # CHROM:POS:REF:ALT from reference vcf
    UNSORTEDQUERY='~{query}'     # CHROM:POS:REF:ALT from query vcf
  
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
  
    PQ="${WORKDIR}/PQ"
    printf -- "${TEMPLATE}" "${PQ}"
    comm -1 -2 "${QUERY}" "${PCA}" > "${PQ}"
    printf -- 'done\n'
  
    if [[ ! -s ${PQ} ]]; then
      error 'No variant in the PCA variants file is mentioned in the query VCF.'
    fi
  
    PR="${WORKDIR}/PR"
    printf -- "${TEMPLATE}" "${PR}"
    comm -1 -2 "${REFERENCE}" "${PCA}" > "${PR}"
    printf -- 'done\n'
  
    if [[ ! -s ${PR} ]]; then
      error 'No variant in the PCA variants file is mentioned in the reference VCF.'
    fi
  
    NPR=$( wc --lines < "${PR}" )
    NPCS=20
    if (( NPR < 2 * NPCS + 1 )); then
  
      # NB: The flashpca program, invoked by PCATasks.PerformPCA, will fail if
      #
      #   min(NVARIANTS, NSAMPLES) < 2 * NPCS + 1
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
  
    if [[ ! -s ${NIX} ]]; then
      exit 0
    fi
  
    (
      sortvariants() {
        local variants="${1}"
        sort                  \
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
        if [[ ${label0} == query ]]; then
          label1=reference
        else
          label1=query
        fi
        if [[ -s ${nixx} ]]; then
          printf -- "${SEPARATOR}"
          if [[ -z "${SEPARATOR}" ]]; then
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

    TEMP_PCA="${WORKDIR}/TEMP_PCA"
    printf -- "${TEMPLATE}" "${TEMP_PCA}"
    comm -2 -3 "${PCA}" "${NIX}" > "${TEMP_PCA}"
    printf -- 'done\n'
  
    WANTED_PCA="${WORKDIR}/WANTED_PCA"
    printf -- "${TEMPLATE}" "${WANTED_PCA}"
    join                                             \
      -t $'\t'                                       \
      "${TEMP_PCA}"                                  \
      <(
        perl -lpe '$_ = qq($_\t$.)' "${UNSORTEDPCA}" \
        | sort                                       \
            --field-separator=$'\t'                  \
            --key=1,1
       )                                             \
      | sort                                         \
        --field-separator=$'\t'                      \
        --key=2,2n                                   \
      | cut --fields=1                               \
      > "${WANTED_PCA}"
    printf -- 'done\n'
  
    mkdir --verbose --parents "$( dirname '~{KEPT_PCA_VARIANTS}' )"
    cp --verbose "${WANTED_PCA}" '~{KEPT_PCA_VARIANTS}'
  
    printf -- '\n\n## WORKDIR:\n'
    find '~{WORKDIR}' -type f | xargs ls -ltr
  
    (
      exec >&2
      printf -- '\n\n### LINE COUNTS:\n'
      find '~{WORKDIR}' -type f | xargs ls -1tr | xargs wc --lines
    )
  
    if ~{if debug then "true" else "false"}; then
      tar -cvzf '~{ARCHIVE}' '~{WORKDIR}'
    fi
  
    printf -- 'FINAL STORAGE UTILIZATION:\n'
    df --human
  
    if ~{if !debug then "true" else "false"}; then
      rm -rf '~{WORKDIR}'
    fi
    >>>
  
  output {
      File?        kept_pca_variants = KEPT_PCA_VARIANTS
      File?        warnings          = WARNINGS
      Array[File?] workfiles         = WORKFILES
      File?        workfiles_tgz     = ARCHIVE
  }
  
  runtime {
      docker: "~{docker_image}"
      disks : "local-disk ~{final_disk_size} HDD"
  }
}  

task BundleAdjustmentModel {
  input {
    Object model_data
    String docker_image
  }

  command {}

  output {
    File manifest = write_json(model_data)
  }

  runtime {
    docker: "~{docker_image}"
  }
}
