version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/Utilities.wdl"
#import v2.2.1
import "https://raw.githubusercontent.com/broadinstitute/warp/f0e6d797fef941c2cfea260a9dd0adcb8effe961/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl"
#import v2.1.12
#import "https://raw.githubusercontent.com/broadinstitute/warp/a4aa63170b09337a3612db6ea22e01b5b332bd54/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl"

#workflow assumes input will be a gvcf file located in a local GCP bucket

workflow ReblockGVCFandPush {
    input {
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"

        # subject, sample id and data location
        String input_gvcf
        String input_gvcf_index
        String rbgvcf_staging_bucket

        # reference genome files
        File ref_dict
        File ref_fasta
        File ref_fasta_index

    }

    #reblock gvcf 
    call ReblockGVCF.ReblockGVCF as Reblock {
        input:
            gvcf = input_gvcf,
            gvcf_index = input_gvcf_index,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            gvcf_file_extension = ".g.vcf.gz",
            cloud_provider = "gcp",
    }


    #Transfer gvcf to staging bucket
    call FileUtils.CopyFilesTask as CopyGVCFToBucket {
    input:
        source_location = Reblock.output_vcf,
        file_types = ["g.vcf.gz"],
        file_match_keys = [],
        file_matchers = [],
        target_location = rbgvcf_staging_bucket,
        flatten = false,
        recursive = true,
        verbose = true,
        docker_image = orchutils_docker_image,
        disk_size = 75,
        gcp_project_id = gcp_project_id,
        workspace_name = workspace_name,
    }

   #Transfer gvcf index to staging bucket
    call FileUtils.CopyFilesTask as CopyGVCFIndexToBucket {
    input:
        source_location = Reblock.output_vcf_index,
        file_types = ["g.vcf.gz.tbi"],
        file_match_keys = [],
        file_matchers = [],
        target_location = rbgvcf_staging_bucket,
        flatten = false,
        recursive = true,
        verbose = true,
        docker_image = orchutils_docker_image,
        disk_size = 75,
        gcp_project_id = gcp_project_id,
        workspace_name = workspace_name,
    }

    output {
        # ReblockGVCF outputs
        File rb_vcf = Reblock.output_vcf
        File rb_vcf_index = Reblock.output_vcf_index
    }
}

