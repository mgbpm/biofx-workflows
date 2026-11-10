version 1.0

import "../../steps/Utilities.wdl"
import "HLATasks.wdl"

workflow T1kGenotypingWorkflow {
    input {
        File    input_sample
        File    input_sample_idx
        File?   hla_intervals
        String? google_project_id
       }

    String format = if (basename(input_sample) == basename(input_sample, ".bam") + ".bam") then "bam"
                    else if (basename(input_sample) == basename(input_sample, ".cram") + ".cram") then "cram"
                    else "unsupported"

    if (format == "unsupported") {
        call Utilities.FailTask as UnsupportedFormatError {
            input:
                error_message = "The input sample file is neither a BAM nor CRAM."
        }
    }

    if (format == "bam" || format == "cram") {
        if (defined(hla_intervals)) {
            call HLATasks.FilterBam as FilterBamToHla {
                input:
                    input_file     = input_sample,
                    input_file_idx = input_sample_idx,
                    intervals      = select_first([hla_intervals]),
                    google_project = google_project_id
            }

            String new_format = "bam"
        }
        
        call HLATasks.RunT1kTask as ImputeHlaAlleles {
            input:
                input_format    = select_first([new_format, format]),
                input_file      = select_first([FilterBamToHla.output_bam, input_sample]),
                input_index     = select_first([FilterBamToHla.output_bai, input_sample_idx]),
                output_basename = sub(basename(input_sample), "\\.(cram|CRAM|bam|BAM)$", ""),
        }

        call HLATasks.PhasingTask as PhaseHaplotypes {
            input:
                input_alleles = ImputeHlaAlleles.allele_tsv
        }
    }

    output {
        File hla_bam         = select_first([FilterBamToHla.output_bam])
        File hla_bai         = select_first([FilterBamToHla.output_bai])
        File haplotypes      = select_first([ImputeHlaAlleles.allele_tsv])
        File haplotype_pairs = select_first([PhaseHaplotypes.haplotype_pairs])
    }
}

