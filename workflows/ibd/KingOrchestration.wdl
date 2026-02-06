version 1.0

import "../../steps/Utilities.wdl"
import "KingTasks.wdl"

workflow KingWorkflow {
    input {
        Array[File]  input_vcfs
        Array[File]  input_vcfs_idx
        File?        input_bed
        String       output_basename
        String       run_type              = "related"
        Int          degree                = 3
        Boolean      missing_to_ref        = false
        Int?         king_mem
        String       bcftools_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        String       king_docker_image     = "uwgac/topmed-master@sha256:0bb7f98d6b9182d4e4a6b82c98c04a244d766707875ddfd8a48005a9f5c5481e"
    }

    if ((run_type != "ibdseg") && (run_type != "related") && (run_type != "kinship") && (run_type != "duplicate")) {
        String RunTypeFail = "Accepted run types include: ibdseg, related, kinship, and duplicate."
    }
    if (length(input_vcfs) != length(input_vcfs_idx)) {
        String IndexCountFail = "The number of index files does not match the number of input VCFs."
    }
    String error_message = select_first([RunTypeFail, IndexCountFail, ""])
    if (error_message != "") {
        call Utilities.FailTask as InputFailure {
            input:
                error_message = error_message
        }
    }

    if (error_message == "") {
        # Merge input VCFs and/or count SNPs and samples
        if (length(input_vcfs) > 1) {
            call KingTasks.MergeVcfsTask as MergeVcfs {
                input:
                    input_vcfs      = input_vcfs,
                    input_vcfs_idx  = input_vcfs_idx,
                    input_bed       = input_bed,
                    missing_to_ref  = missing_to_ref,
                    output_basename = output_basename,
                    docker_image    = bcftools_docker_image
            }
        }
        if (length(input_vcfs) == 1) {
            File dataset_vcf     = input_vcfs[0]
            File dataset_vcf_idx = input_vcfs_idx[0]

            call KingTasks.FilterVcfTask as PrepSingleVcf {
                input:
                    input_vcf     = dataset_vcf,
                    input_vcf_idx = dataset_vcf_idx,
                    input_bed     = input_bed,
                    docker_image  = bcftools_docker_image
            }
        }
        File input_vcf_ = select_first([
            MergeVcfs.merged_vcf,
            PrepSingleVcf.output_vcf_gz,
            dataset_vcf
        ])

        # Convert all-sample VCF to PLINK files
        call KingTasks.Vcf2BedTask as Vcf2Bed {
            input:
                input_vcf       = input_vcf_,
                output_basename = output_basename,
                docker_image    = king_docker_image
        }

        # Estimate memory needed for KING
        Int num_samples  = select_first([MergeVcfs.num_samples, PrepSingleVcf.num_samples])
        Int num_snps     = select_first([MergeVcfs.num_snps, PrepSingleVcf.num_snps])
        Int est_king_mem = ceil((num_samples * (num_snps/4)) / 1000000000)

        Int king_mem_ = if est_king_mem < 2 then 4 else est_king_mem + 4
        Int final_king_mem = select_first([king_mem, king_mem_])

        # Run KING with specified flag
        call KingTasks.RunKingTask {
            input:
                bed_file        = Vcf2Bed.bed_file,
                bim_file        = Vcf2Bed.bim_file,
                fam_file        = Vcf2Bed.fam_file,
                flag            = run_type,
                degree          = degree,
                output_basename = output_basename,
                docker_image    = king_docker_image,
                mem_size        = final_king_mem
        }
    }

    output {
        File? seg_output  = RunKingTask.seg_output
        File? con_output  = RunKingTask.con_output
        File? kin_output  = RunKingTask.kin_output
        File? kin0_output = RunKingTask.kin0_output
    }
}