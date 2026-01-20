version 1.0

import "../../steps/Utilities.wdl"
import "KingTasks.wdl"

workflow KingWorkflow {
    input {
        Array[File] input_vcfs
        Array[File]? input_vcfs_idx
        String output_basename
        File? input_bed
        String run_type = "related"
        Int degree = 3
        Int? king_mem
        String bcftools_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        String king_docker_image = "uwgac/topmed-master@sha256:0bb7f98d6b9182d4e4a6b82c98c04a244d766707875ddfd8a48005a9f5c5481e"
    }

    if ((run_type != "ibdseg") && (run_type != "related") && (run_type != "kinship") && (run_type != "duplicate")) {
        call Utilities.FailTask as RunTypeFail {
            input:
                error_message = "Incorrect input for run_type input parameter."
        }
    }

    # Merge all input VCFs
    if (length(input_vcfs) > 1) {
        if (!defined(input_vcfs_idx)) {
            call Utilities.FailTask as IndexFail {
                input:
                    error_message = "Include VCF index files for merging."
            }
        }
        if (length(input_vcfs) != length(select_first([input_vcfs_idx]))) {
            call Utilities.FailTask as NumIndexFail {
                input:
                    error_message = "The number of index files does not match the number of input VCFs."
            }
        }

        call KingTasks.MergeVcfsTask as MergeVcfs {
            input:
                input_vcfs = input_vcfs,
                input_vcfs_idx = select_first([input_vcfs_idx]),
                input_bed = input_bed,
                output_basename = output_basename,
                docker_image = bcftools_docker_image
        }
    }
    # Skip merging if there is only one input VCF
    if (length(input_vcfs) == 1) {
        File dataset_vcf = input_vcfs[0]
        File dataset_vcf_idx = select_first([input_vcfs_idx])[0]

        # Even if no BED file is provided, this will count SNPs and samples
        call KingTasks.FilterVcfTask as FilterSingleVcf {
            input:
                input_vcf = dataset_vcf,
                input_bed = input_bed,
                docker_image = bcftools_docker_image
        }
    }
    File input_vcf_ = select_first([
        MergeVcfs.merged_vcf,
        FilterSingleVcf.output_vcf_gz,
        dataset_vcf
    ])

    # Convert merged VCF to PLINK bed
    call KingTasks.Vcf2BedTask as Vcf2Bed {
        input:
            input_vcf = input_vcf_,
            output_basename = output_basename,
            docker_image = king_docker_image
    }

    # Estimate memory needed for KING
    Int num_samples = select_first([MergeVcfs.num_samples, FilterSingleVcf.num_samples])
    Int num_snps = select_first([MergeVcfs.num_snps, FilterSingleVcf.num_snps])
    Int est_king_mem = ceil((num_samples * (num_snps/4)) / 1000000000)

    if (est_king_mem < 2) {
        Int king_mem__ = 4
    }
    if (est_king_mem > 2) {
        Int king_mem_ = est_king_mem + 4
    }
    Int final_king_mem = select_first([king_mem, king_mem__, king_mem_])

    # Run KING with specified flag
    call KingTasks.RunKingTask {
        input:
            bed_file = Vcf2Bed.bed_file,
            bim_file = Vcf2Bed.bim_file,
            fam_file = Vcf2Bed.fam_file,
            flag = run_type,
            degree = degree,
            output_basename = output_basename,
            docker_image = king_docker_image,
            mem_size = final_king_mem
    }

    output {
        File? seg_output = RunKingTask.seg_output
        File? con_output = RunKingTask.con_output
        File? kin_output = RunKingTask.kin_output
        File? kin0_output = RunKingTask.kin0_output
    }
}