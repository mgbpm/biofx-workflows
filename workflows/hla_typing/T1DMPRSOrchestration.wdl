version 1.0

import "./T1KGenotyping.wdl"
import "./HLATasks.wdl"
import "./RawScoring.wdl"
import "../lowpassimputation/Glimpse2Imputation.wdl"

workflow T1dmPrsWorkflow {
    input {
        File    input_cram
        File    input_crai
        File?   input_vcf
        File?   input_pgen
        File?   input_pvar
        File?   input_psam
        File    t1dm_weights_file
        File?   hla_intervals
        String? google_project_id
        File    ref_fasta                    = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File    ref_fai                      = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        File    ref_dict                     = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
        File    glimpse_reference_chunks     = "gs://lmm-reference-data/prsmix/glimpse_reference_chunks/hgdp_1kg_ref_chunks_withChrX/list/hgdp_1kg_ref_chunks_withChrX"
        File    glimpse_ref_vcf              = "gs://lmm-reference-data/prsmix/reference/NA12878_437224286_v1.1.imputed_VEP_split_af_sorted.vcf.gz"
        String  glimpse_docker_image         = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
        String  glimpse_extract_docker_image = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
    }

    String  sample_basename = sub(basename(input_cram), "\\.(cram|CRAM|bam|BAM)$", "")
    Boolean has_pgen_input  = defined(input_pgen) || defined(input_pvar) || defined(input_psam)
    Boolean has_vcf_input   = defined(input_vcf)

    if (!has_vcf_input && !has_pgen_input) {
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

    File scoring_vcf = select_first([input_vcf, RunGlimpse.imputed_afFiltered_vcf])

    call RawScoring.RawScoreWorkflow as ScorePrs {
        input:
            input_vcf           = scoring_vcf,
            input_pgen          = input_pgen,
            input_pvar          = input_pvar,
            input_psam          = input_psam,
            weights             = t1dm_weights_file
    }

    call T1KGenotyping.T1kGenotypingWorkflow as GetHlaHaplotypes {
        input:
            input_sample      = input_cram,
            input_sample_idx  = input_crai,
            hla_intervals     = select_first([hla_intervals]),
            google_project_id = select_first([google_project_id])
    }

    call HLATasks.ScoreHaplotypesTask as ScoreHlaHaplotypes {
        input:
            haplotype_pairs = GetHlaHaplotypes.haplotype_pairs
    }

    call MixScoresTask as MixHlaAndPrs {
        input:
            hla_score_file = ScoreHlaHaplotypes.scored_haplotypes,
            prs_score_file = ScorePrs.prs_score,
            sample_id      = sample_basename
    }

    output {
        File hla_haplotype_pairs = GetHlaHaplotypes.haplotype_pairs
        File hla_score           = ScoreHlaHaplotypes.scored_haplotypes
        File vcf                 = scoring_vcf
        File prs_score           = ScorePrs.prs_score
        File final_score         = MixHlaAndPrs.mixed_score
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
        weighted_prs_score=$(awk -v x=8.86562 -v y=${prs_score} 'BEGIN {print x*y}')

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