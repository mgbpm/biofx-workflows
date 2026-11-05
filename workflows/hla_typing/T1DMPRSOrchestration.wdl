version 1.0

import "./T1KGenotyping.wdl"
import "./HLAPhasing.wdl"
import "./RawScoring.wdl"
import "../../steps/VCFUtils.wdl"
import "../lowpassimputation/Glimpse2Imputation.wdl"

workflow T1dmPrsWorkflow {
    input {
        File    input_cram
        File    input_crai
        File?   input_vcf
        File    t1dm_weights_file
        String  vcf_generation               = "glimpse"
        File?   gatk_roi
        File?   hla_intervals
        String? google_project_id
        File    ref_fasta                    = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File    ref_fai                      = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        File    ref_dict                     = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
        File    glimpse_reference_chunks     = "gs://lmm-reference-data/prsmix/glimpse_reference_chunks/hgdp_1kg_ref_chunks_withChrX/list/hgdp_1kg_ref_chunks_withChrX"
        File    glimpse_ref_vcf              = "gs://lmm-reference-data/prsmix/reference/NA12878_437224286_v1.1.imputed_VEP_split_af_sorted.vcf.gz"
        String  t1k_docker_image             = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/t1k:0.0.1"
        String  prs_docker_image             = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:dev-hla"
        String  bcftools_docker_image        = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        String  plink2_docker_image          = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        String  glimpse_docker_image         = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
        String  glimpse_extract_docker_image = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
        String  gatk_docker_image             = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    }

    String sample_basename = sub(basename(input_cram), "\\.(cram|CRAM|bam|BAM)$", "")

    call T1KGenotyping.T1kGenotypingWorkflow as GetHlaHaplotypes {
        input:
            input_sample      = input_cram,
            input_sample_idx  = input_crai,
            hla_intervals     = select_first([hla_intervals]),
            google_project_id = select_first([google_project_id]),
            ref_fasta         = ref_fasta,
            ref_fai           = ref_fai,
            ref_dict          = ref_dict,
            t1k_docker_image  = t1k_docker_image
    }

    call HLATasks.ScoreHaplotypesTask as ScoreHlaHaplotypes {
        input:
            haplotype_pairs = GetHlaHaplotypes.haplotype_pairs
    }

    if (!defined(input_vcf)) {
        if (vcf_generation == "gatk") {
            call HaplotypeCallerTask as RunGATK {
                input:
                    input_file      = input_cram,
                    input_file_idx  = input_crai,
                    ref_fasta       = ref_fasta,
                    ref_fai         = ref_fai,
                    ref_dict        = ref_dict,
                    roi_bed         = select_first([gatk_roi]),
                    docker_image    = gatk_docker_image,
                    gatk_path       = "/gatk/gatk",
                    output_filename = sample_basename + ".vcf.gz"
            }

            call VCFUtils.SortVCFTask as SortGatkVcf {
                input:
                    input_vcf    = RunGATK.output_vcf,
                    docker_image = bcftools_docker_image

            }
        }

        if (vcf_generation == "glimpse") {
            call Glimpse2Imputation.Glimpse2Imputation as RunGlimpse {
                input:
                    reference_chunks                              = glimpse_reference_chunks,
                    crams                                         = [input_cram],
                    cram_indices                                  = [input_crai],
                    sample_ids                                    = [sample_basename],
                    fasta                                         = ref_fasta,
                    fasta_index                                   = ref_fai,
                    ref_dict                                      = ref_dict,
                    output_basename                               = sample_basename,
                    af_cutoff                                     = ">=0.0001",
                    gnomadAF_ref_vcf                              = glimpse_ref_vcf,
                    docker                                        = glimpse_docker_image,
                    docker_extract_num_sites_from_reference_chunk = glimpse_extract_docker_image
            }
        }
    }

    call RawScoring.RawScoreWorkflow as ScorePrs {
        input:
            input_vcf           = select_first([input_vcf, RunGlimpse.imputed_afFiltered_vcf, SortGatkVcf.output_vcf_gz]),
            weights             = t1dm_weights_file,
            plink2_docker_image = plink2_docker_image
    }

    call MixScoresTask as MixHlaAndPrs {
        input:
            hla_score_file = ScoreHlaHaplotypes.scored_haplotypes,
            prs_score_file = ScorePrs.prs_score,
            sample_id      = sample_basename
    }

    output {
        File hla_genotypes = ImputeHlaGenotype.allele_tsv
        File hla_score     = PhaseHlaHaplotypes.scored_haplotypes
        File vcf           = select_first([input_vcf, RunGlimpse.imputed_afFiltered_vcf, SortGatkVcf.output_vcf_gz])
        File prs_score     = ScorePrs.prs_score
        File final_score   = MixHlaAndPrs.mixed_score
    }
}

task HaplotypeCallerTask {
    input {
        File   input_file
        File   input_file_idx
        File?  roi_bed
        String output_filename
        File   ref_dict
        File   ref_fasta
        File   ref_fai
        String gatk_path
        String docker_image
        Int    addldisk    = 30
        Int    preemptible = 1
        Int    mem         = 24
        Int    bootdisk    = 15
    }

    Int ref_size        = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB") + size(ref_dict, "GB"))
    Int input_size      = ceil(size(input_file, "GB") + size(input_file_idx, "GB") + size(roi_bed, "GB"))
    Int final_disk_size = ref_size + input_size + addldisk

    command <<<
        set -euxo pipefail

        BED_FILE=""
        [ ! -z "~{roi_bed}" ] && BED_FILE="--intervals ~{roi_bed}"

        ~{gatk_path} --java-options "-Xmx~{mem-1}G"       \
            HaplotypeCaller                               \
            --reference ~{ref_fasta}                      \
            --input ~{input_file}                         \
            ${BED_FILE}                                   \
            --output ~{output_filename}
    >>>

    runtime {
        docker: "~{docker_image}"
        memory: "~{mem} GB"
        disks: "local-disk ~{final_disk_size} SSD"
        preemptible: preemptible
        bootDiskSizeGb: bootdisk
    }

    output {
        File output_vcf     = "~{output_filename}"
        File output_vcf_idx = "~{output_filename}.tbi"
    }
}

task MixScoresTask {
    input {
        File   hla_score_file
        File   prs_score_file
        String sample_id
        String docker_image = "ubuntu:latest"
        Int    addldisk = 10
        Int    mem_size = 2
    }

    Int input_size = ceil(size(hla_score_file, "GB") + size(prs_score_file, "GB"))
    Int disk_size  = input_size + addldisk

    command <<<
        set -euxo pipefail

        hla_score=$(grep "~{sample_id}" "~{hla_score_file}" | cut -f14)
        weighted_hla_score=$(awk -v x=0.79003 -v y=${hla_score} 'BEGIN {print x*y}')

        prs_score=$(grep "~{sample_id}" "~{prs_score_file}" | cut -f4)
        weighted_t1dm_score=$(awk -v x=08.86562 -v y=${hla_score} 'BEGIN {print x*y}')

        merged_score=$(awk -v x=${weighted_hla_score} -v y=${weighted_prs_score} 'BEGIN {print x+y}')
        final_score=$(awk -v x=20 -v y=${merged_score} 'BEGIN {print x+y}')

        printf -- "${final_score}" > "~{sample_id}.score.tsv"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_size} GB"
    }

    output {
        File mixed_score = "~{sample_id}.score.tsv"
    }
}