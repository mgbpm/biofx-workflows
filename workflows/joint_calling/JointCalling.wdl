version 1.0

#import "../../steps/FileUtils.wdl"
#import "../../steps/Utilities.wdl"
#import "../../steps/GenerateSampleMap.wdl"
#import "../../steps/JointGenotyping.wdl"
import "https://raw.githubusercontent.com/broadinstitute/warp/32bcffefba4f0aa9671a321dbd076fbc5b616dbc/pipelines/broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl" as JointGenotyping
import "https://raw.githubusercontent.com/gatk-workflows/utility-wdls/main/generate-sample-map.wdl" as GenerateSampleMap


workflow JointCalling {
    input {
        #define GenerateSampleMap inputs
        Array[String] sample_names
        Array[String] file_paths
        String sample_map_name

        #define JointGenotyping inputs
        File unpadded_intervals_file

        #String callset_name
        #File sample_name_map

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File dbsnp_vcf
        File dbsnp_vcf_index

        Int small_disk
        Int medium_disk
        Int large_disk
        Int huge_disk

        Array[String] snp_recalibration_tranche_values
        Array[String] snp_recalibration_annotation_values
        Array[String] indel_recalibration_tranche_values
        Array[String] indel_recalibration_annotation_values

        File haplotype_database

        File eval_interval_list
        File hapmap_resource_vcf
        File hapmap_resource_vcf_index
        File omni_resource_vcf
        File omni_resource_vcf_index
        File one_thousand_genomes_resource_vcf
        File one_thousand_genomes_resource_vcf_index
        File mills_resource_vcf
        File mills_resource_vcf_index
        File axiomPoly_resource_vcf
        File axiomPoly_resource_vcf_index
        File dbsnp_resource_vcf = dbsnp_vcf
        File dbsnp_resource_vcf_index = dbsnp_vcf_index

        # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
        # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
        Float excess_het_threshold = 54.69
        Float? vqsr_snp_filter_level
        Float? vqsr_indel_filter_level
        Int? snp_vqsr_downsampleFactor

        Int? top_level_scatter_count
        Boolean? gather_vcfs
        Int snps_variant_recalibration_threshold = 500000
        Boolean rename_gvcf_samples = true
        Float unbounded_scatter_count_scale_factor = 0.15
        Int gnarly_scatter_count = 10
        Boolean use_gnarly_genotyper = false
        Boolean use_allele_specific_annotations = true
        Boolean cross_check_fingerprints = true
        Boolean scatter_cross_check_fingerprints = false
    }

    call GenerateSampleMap.GenerateSampleMap {
        input:
            file_paths = file_paths,
            sample_map_name = sample_map_name,
            sample_names = sample_names
    }

    File sample_map = GenerateSampleMap.sample_map

    call JointGenotyping.JointGenotyping {
        input:
            unpadded_intervals_file = unpadded_intervals_file,
            callset_name = sample_map_name,
            sample_name_map = GenerateSampleMap.sample_map,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            dbsnp_vcf = dbsnp_vcf,
            dbsnp_vcf_index = dbsnp_vcf_index,
            small_disk = small_disk,
            medium_disk = medium_disk,
            large_disk = large_disk,
            huge_disk = huge_disk,
            snp_recalibration_tranche_values = snp_recalibration_tranche_values,
            snp_recalibration_annotation_values = snp_recalibration_annotation_values,
            indel_recalibration_tranche_values = indel_recalibration_tranche_values,
            indel_recalibration_annotation_values = indel_recalibration_annotation_values,
            haplotype_database = haplotype_database,
            eval_interval_list = eval_interval_list,
            hapmap_resource_vcf = hapmap_resource_vcf,
            hapmap_resource_vcf_index = hapmap_resource_vcf,
            omni_resource_vcf = omni_resource_vcf,
            omni_resource_vcf_index = omni_resource_vcf_index,
            one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
            one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
            mills_resource_vcf = mills_resource_vcf,
            mills_resource_vcf_index = mills_resource_vcf_index,
            axiomPoly_resource_vcf = axiomPoly_resource_vcf,
            axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
            dbsnp_resource_vcf = dbsnp_vcf,
            dbsnp_resource_vcf_index = dbsnp_vcf_index,

            # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
            # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
            excess_het_threshold = 54.69,
            
            vqsr_snp_filter_level = vqsr_snp_filter_level,
            vqsr_indel_filter_level = vqsr_indel_filter_level,
            snp_vqsr_downsampleFactor = snp_vqsr_downsampleFactor,

            top_level_scatter_count = top_level_scatter_count,
            gather_vcfs = gather_vcfs,
            
            snps_variant_recalibration_threshold = 500000,
            rename_gvcf_samples = true,
            unbounded_scatter_count_scale_factor = 0.15,
            gnarly_scatter_count = 10,
            use_gnarly_genotyper = false,
            use_allele_specific_annotations = true,
            cross_check_fingerprints = true,
            scatter_cross_check_fingerprints = false
    }
    output {
        # Metrics from either the small or large callset
        File detail_metrics_file = JointGenotyping.detail_metrics_file
        File summary_metrics_file = JointGenotyping.summary_metrics_file

        # Outputs from the small callset path through the wdl.
        Array[File] output_vcfs = select_all(JointGenotyping.output_vcfs)
        Array[File] output_vcf_indices = select_all(JointGenotyping.output_vcf_indices)

        # Output the interval list generated/used by this run workflow.
        Array[File] output_intervals = JointGenotyping.output_intervals

        # Output the metrics from crosschecking fingerprints.
        File? crosscheck_fingerprint_check = JointGenotyping.crosscheck_fingerprint_check
    }
}