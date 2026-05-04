version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/HaplotypeCallerGvcfGATK4.wdl"
import "../pgxrisk/PGxWorkflow.wdl"
import "../pgxrisk/RiskAllelesWorkflow.wdl"
import "../../steps/DepthOfCoverage.wdl"
import "../../steps/VCFUtils.wdl"
import "../../steps/GDHIngestAndFilter.wdl"
import "../../steps/QCEval.wdl"
import "../../steps/IgvReport.wdl"
import "../../steps/Utilities.wdl"

workflow BgwgsWorkflow {
    input {
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
        # bcftools docker image
        String bcftools_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        # subject, sample id and data location
        String subject_id
        String sample_id
        String sample_data_location
        Boolean fetch_cram = true
        Array[String] fetch_cram_filter_keys = [subject_id, sample_id]
        Array[FileMatcher]? fetch_cram_file_matchers
        Boolean fetch_bam = true
        Array[String] fetch_bam_filter_keys = [subject_id, sample_id]
        Array[FileMatcher]? fetch_bam_file_matchers
        Array[String] fetch_vcf_filter_keys = [subject_id, sample_id]
        Array[FileMatcher]? fetch_vcf_file_matchers
        Boolean fetch_files_verbose = false
        Boolean do_variant_calling = true
        Boolean do_coverage = true
        Boolean do_pgx = true
        Boolean do_risk_alleles = true
        Boolean do_gnomad = true
        Int fetch_disk_size = 75
        # reference genome files
        String reference_build = "GRCh38"
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        File? dbsnp_vcf
        File? dbsnp_vcf_index
        # variant calling inputs
        File? scattered_calling_intervals_list
        # depth of coverage inputs
        File cov_roi_bed = "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed"
        Array[RoiAndRefGeneFilePair?] cov_roi_genes = [
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi0.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene0.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene0.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi1.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene1.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene1.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi10.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene10.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene10.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi11.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene11.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene11.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi12.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene12.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene12.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi13.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene13.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene13.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi14.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene14.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene14.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi15.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene15.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene15.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi16.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene16.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene16.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi17.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene17.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene17.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi18.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene18.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene18.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi19.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene19.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene19.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi2.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene2.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene2.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi20.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene20.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene20.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi21.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene21.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene21.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi22.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene22.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene22.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi24.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene24.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene24.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi25.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene25.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene25.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi26.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene26.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene26.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi27.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene27.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene27.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi29.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene29.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene29.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi3.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene3.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene3.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi30.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene30.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene30.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi31.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene31.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene31.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi32.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene32.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene32.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi34.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene34.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene34.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi36.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene36.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene36.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi4.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene4.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene4.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi5.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene5.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene5.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi6.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene6.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene6.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi7.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene7.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene7.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi8.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene8.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene8.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi9.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene9.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene9.txt.idx"
            },
            {
                "roi_bed": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/roi97.bed",
                "ref_gene": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene97.txt",
                "ref_gene_idx": "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed_splitoverlap/grp/refgene97.txt.idx"
            }
        ]
        File cov_gene_names = "gs://lmm-reference-data/roi/HGNC_genenames_05272022.txt"
        String cov_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/coverage:20230630"
        String gatk3_docker_image = "broadinstitute/gatk3:3.7-0"
        # pgx inputs
        String pgx_test_code = "lmPGX-pnlD_L"
        String pgx_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20241008"
        File pgx_workflow_fileset = "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_20241004.tar"
        File pgx_roi_bed = "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_genotyping.bed"
        # risk alleles inputs
        String risk_alleles_test_code = "lmRISK-pnlB_L"
        String risk_alleles_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20240129"
        File risk_alleles_workflow_fileset = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar"
        File risk_alleles_roi_bed = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed"
        # qceval inputs
        String qceval_project_type = "BGE_DRAGEN_TP_BINNING"
        String qceval_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20250923"
        File? thresholds = "gs://lmm-reference-data/annotation/pmeval/thresholds_20250912.tsv"
        File? difficult_to_map_regions = "gs://lmm-reference-data/annotation/pmeval/difficult_to_map_regions_20250912.tgz"
        # gnomad annotation inputs
        File gnomad_coverage_file = "gs://lmm-reference-data/annotation/gnomad/genomes.r3.0.1.coverage_targetROI-filtered.dedup.txt.gz"
        File gnomad_coverage_file_idx = "gs://lmm-reference-data/annotation/gnomad/genomes.r3.0.1.coverage_targetROI-filtered.dedup.txt.gz.tbi"
        Array[String] gnomad_headers = [ "##INFO=<ID=DP_gnomadG,Number=1,Type=Float,Description=\"Read depth of GnomAD Genome\">" ]
        String gnomad_column_list = "CHROM,POS,INFO/DP_gnomadG"
        # GDH ingest and filter
        String gdh_institution = "MGBPM"
        String gdh_project = "Clinical"
        String vcf_file_stage_name = "biofx_pipelines"
        String vcf_file_stage_gspath = "gs://gdh-external-stage/biofx_pipelines_nonprod"
        String filter_name_or_code
        String? pipeline_run_id
        # Reporting steps
        String igvreport_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511"
        String gdh_parser_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/gdhoutputparser:dev"
        File portable_db_file = "gs://lmm-reference-data/annotation/gil_lmm/gene_info.db"
        Array[File] igv_track_files = [ "gs://lmm-reference-data/annotation/ucsc/hg38/refGene_20231019.txt.gz" ]
        Array[File] igv_track_index_files = [ "gs://lmm-reference-data/annotation/ucsc/hg38/refGene_20231019.txt.gz.tbi" ]
    }

    # Prefer a BAM file to avoid conversion
    if (fetch_bam) {
        call FileUtils.FetchFilesTask as FetchBam {
            input:
                data_location = sample_data_location,
                file_types = if defined(fetch_bam_file_matchers) then [] else ["bam"],
                recursive = false,
                file_match_keys = if defined(fetch_bam_file_matchers) then [] else fetch_bam_filter_keys,
                file_matchers = fetch_bam_file_matchers,
                verbose = fetch_files_verbose,
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                disk_size = fetch_disk_size
        }
    }

    # no BAM, fallback to CRAM
    if (fetch_cram && (!defined(FetchBam.bam) || !defined(FetchBam.bai))) {
        call FileUtils.FetchFilesTask as FetchCram {
            input:
                data_location = sample_data_location,
                file_types = if defined(fetch_cram_file_matchers) then [] else ["cram"],
                recursive = false,
                file_match_keys = if defined(fetch_cram_file_matchers) then [] else fetch_cram_filter_keys,
                file_matchers = fetch_cram_file_matchers,
                verbose = fetch_files_verbose,
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                disk_size = fetch_disk_size
        }
    }

    # Validate that we got workable files from the data location
    if (!defined(FetchBam.bam) && !defined(FetchCram.cram)) {
        call Utilities.FailTask as MissingBamOrCramFailure {
            input:
                error_message = "BAM or CRAM file not found in ~{sample_data_location}"
        }
    }
    if (defined(FetchCram.cram) && !defined(FetchCram.crai)) {
        call Utilities.FailTask as MissingCraiFailure {
            input:
                error_message = "Index file for CRAM " + basename(select_first([FetchCram.cram])) + " not found"
        }
    }
    if (!defined(FetchCram.cram) && defined(FetchBam.bam) && !defined(FetchBam.bai)) {
        call Utilities.FailTask as MissingBaiFailure {
            input:
                error_message = "Index file for BAM " + basename(select_first([FetchCram.bam])) + " not found"
        }
    }

    # Convert CRAM to BAM
    if (defined(FetchCram.cram) && defined(FetchCram.crai)) {
        call HaplotypeCallerGvcfGATK4.CramToBamTask {
            input:
                input_cram = select_first([FetchCram.cram]),
                sample_name = basename(select_first([FetchCram.cram]), ".cram"),
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710",
                samtools_path = "samtools"
        }
    }

    # final inputs - either CRAM or BAM
    File sample_bam = select_first([FetchBam.bam, CramToBamTask.output_bam])
    File sample_bai = select_first([FetchBam.bai, CramToBamTask.output_bai])

    # Run haplotype caller
    if (do_variant_calling) {
        call HaplotypeCallerGvcfGATK4.HaplotypeCallerGvcf_GATK4 {
            input:
                input_bam = sample_bam,
                input_bam_index = sample_bai,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                scattered_calling_intervals_list = select_first([scattered_calling_intervals_list]),
                make_gvcf = false
        }
    }
    if (!do_variant_calling) {
        # Fetch the VCF files from the storage location
        call FileUtils.FetchFilesTask as FetchVcf {
            input:
                data_location = sample_data_location,
                file_types = if defined(fetch_vcf_file_matchers) then [] else ["vcf"],
                recursive = false,
                file_match_keys = if defined(fetch_vcf_file_matchers) then [] else fetch_vcf_filter_keys,
                file_matchers = fetch_vcf_file_matchers,
                verbose = fetch_files_verbose,
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name
        }
        if (!defined(FetchVcf.vcf)) {
            call Utilities.FailTask as MissingVcfFailure {
                input:
                    error_message = "VCF file not found in ~{sample_data_location}"
            }
        }
    }

    File sample_vcf = select_first([HaplotypeCallerGvcf_GATK4.output_vcf, FetchVcf.vcf])

    # Run depth of coverage
    if (do_coverage) {
        call DepthOfCoverage.DepthOfCoverageWorkflow {
            input:
                run_wgs = true,
                output_basename = subject_id + "_" + sample_id + ".cov",
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bam = sample_bam,
                bai = sample_bai,
                roi_all_bed = cov_roi_bed,
                roi_genes = cov_roi_genes,
                gene_names = cov_gene_names,
                cov_docker_image = cov_docker_image,
                gatk_docker_image = gatk3_docker_image
        }
    }

    # Run PGx
    if (do_pgx) {
        call PGxWorkflow.PGxWorkflow as PGxWorkflowAlias {
            input:
                input_cram = sample_bam,
                input_crai = sample_bai,
                sample_id = sample_id,
                accession_id = subject_id,
                test_code = pgx_test_code,
                reference_fasta = ref_fasta,
                reference_fasta_fai = ref_fasta_index,
                reference_dict = ref_dict,
                roi_bed = pgx_roi_bed,
                dbsnp = select_first([dbsnp_vcf]),
                dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
                workflow_fileset = pgx_workflow_fileset,
                mgbpmbiofx_docker_image = pgx_docker_image
        }
    }

    # Run Risk
    if (do_risk_alleles) {
        call RiskAllelesWorkflow.RiskAllelesWorkflow as RiskAllelesWorkflowAlias {
            input:
                input_cram = sample_bam,
                input_crai = sample_bai,
                sample_id = sample_id,
                accession_id = subject_id,
                test_code = risk_alleles_test_code,
                reference_fasta = ref_fasta,
                reference_fasta_fai = ref_fasta_index,
                reference_dict = ref_dict,
                roi_bed = risk_alleles_roi_bed,
                dbsnp = select_first([dbsnp_vcf]),
                dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
                workflow_fileset = risk_alleles_workflow_fileset,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                mgbpmbiofx_docker_image = risk_alleles_docker_image
        }
    }


    # The assignments at the end of the conditional block below ensure
    # that `ref_fasta` and `ref_fasta_index` are passed to
    # `QCEval.QCEvalTask` only when needed, avoiding unnecessary
    # localization--especially of `ref_fasta`, which is typically
    # large.  (Because these inputs have type `File`, they are
    # mandatory and cannot be omitted.  If passed directly to the
    # task, they would be localized unconditionally.  Assigning them
    # via a conditional block adds them to the task's inputs only when
    # needed.  In contrast, the other inputs specific to the
    # "BGE_DRAGEN_TP_BINNING" case--`thresholds` and
    # `difficult_to_map_regions`--have type `File?`, so they can
    # simply be left out when not needed; localization is not an issue
    # in this case.)

    if (qceval_project_type == "BGE_DRAGEN_TP_BINNING") {
      if (!(defined(thresholds) && defined(difficult_to_map_regions))) {
        call Utilities.FailTask as MissingQcevalBgeDragenTpBinningInputs {
            input:
                error_message = "Some inputs required for QCEval with qceval_project_type = \\\"~{qceval_project_type}\\\" are missing"
        }
      }
      File maybe_reference_fasta = ref_fasta
      File maybe_reference_fasta_fai = ref_fasta_index
    }

    # Annotate target VCF with QCEval INFO field
    call QCEval.QCEvalTask {
        input:
            input_vcf = sample_vcf,
            project_type = qceval_project_type,
            output_basename = subject_id + "_" + sample_id + ".qceval",
            reference_fasta = maybe_reference_fasta,
            reference_fasta_fai = maybe_reference_fasta_fai,
            thresholds_tsv = thresholds,
            regions_tgz = difficult_to_map_regions,
            docker_image = qceval_docker_image
    }

    # Annotated target VCF with Gnomad coverage
    if (do_gnomad) {
        call VCFUtils.AnnotateVCFTask as AnnotateGnomadTask {
            input:
                input_vcf = QCEvalTask.output_vcf_gz,
                output_basename = subject_id + "_" + sample_id + ".gnomad",
                annotations_file = gnomad_coverage_file,
                annotations_idx_file = gnomad_coverage_file_idx,
                headers_file = write_lines(gnomad_headers),
                column_list = gnomad_column_list,
                docker_image = "~{bcftools_docker_image}"
        }
    }

    # Load sample data to GDH and run filtration
    call GDHIngestAndFilter.GDHIngestAndFilterTask {
        input:
            subject_id = subject_id,
            sample_id = sample_id,
            gdh_institution = gdh_institution,
            gdh_project = gdh_project,
            vcf_file = select_first([AnnotateGnomadTask.output_vcf_gz, QCEvalTask.output_vcf_gz]),
            vcf_file_stage_name = vcf_file_stage_name,
            vcf_file_stage_gspath = vcf_file_stage_gspath,
            reference_build = reference_build,
            filter_name_or_code = filter_name_or_code,
            pipeline_run_id = pipeline_run_id,
            timeout_minutes = 90,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }

    call GDHIngestAndFilter.GDHOutputParserTask {
        input:
            gdh_output_file = GDHIngestAndFilterTask.matching_variants,
            reference_build = reference_build,
            oms_query = "Y",
            portable_db_file = portable_db_file,
            filter_name_or_code = filter_name_or_code,
            report_basename = "~{subject_id}_~{sample_id}_~{filter_name_or_code}",
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            gdh_parser_image = gdh_parser_image
    }
    if (defined(GDHOutputParserTask.parsed_report)) {
        call IgvReport.IgvReportFromParsedFASTOutputTask {
            input:
                bam_cram = sample_bam,
                bai_crai = sample_bai,
                parsed_fast_output = select_first([GDHOutputParserTask.parsed_report]),
                output_basename = "~{subject_id}_~{sample_id}_~{filter_name_or_code}.igvreport",
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                track_files = igv_track_files,
                track_index_files = igv_track_index_files,
                docker_image = igvreport_docker_image
        }
    }
    call BgwgsGDHSummaryTask {
        input:
            subject_id = subject_id,
            sample_id = sample_id,
            filter_name_or_code = filter_name_or_code
    }

    output {
        # coverage outputs
        File? cov_wgs_sample_summary = DepthOfCoverageWorkflow.wgs_sample_summary
        File? cov_wgs_sample_statistics = DepthOfCoverageWorkflow.wgs_sample_statistics
        File? cov_roi_sample_interval_summary = DepthOfCoverageWorkflow.roi_sample_interval_summary
        File? cov_roi_sample_interval_statistics = DepthOfCoverageWorkflow.roi_sample_interval_statistics
        File? cov_roi_sample_statistics = DepthOfCoverageWorkflow.roi_sample_statistics
        File? cov_roi_sample_summary = DepthOfCoverageWorkflow.roi_sample_summary
        File? cov_roi_sample_cumulative_coverage_counts = DepthOfCoverageWorkflow.roi_sample_cumulative_coverage_counts
        File? cov_roi_sample_cumulative_coverage_proportions = DepthOfCoverageWorkflow.roi_sample_cumulative_coverage_proportions
        File? cov_mt_summary = DepthOfCoverageWorkflow.mt_summary
        File? cov_gene_summary = DepthOfCoverageWorkflow.gene_summary
        File? cov_gene_summary_unknown = DepthOfCoverageWorkflow.gene_summary_unknown
        File? cov_gene_summary_entrez = DepthOfCoverageWorkflow.gene_summary_entrez
        # haplotype caller output
        File? vcf = HaplotypeCallerGvcf_GATK4.output_vcf
        # pgx output
        File? pgx_summary_report = PGxWorkflowAlias.summary_report
        File? pgx_details_report = PGxWorkflowAlias.details_report
        File? pgx_genotype_xlsx = PGxWorkflowAlias.genotype_xlsx
        File? pgx_genotype_txt = PGxWorkflowAlias.genotype_txt
        # risk alleles output
        File? risk_alleles_report = RiskAllelesWorkflowAlias.risk_report
        File? risk_alleles_genotype_xlsx = RiskAllelesWorkflowAlias.genotype_xlsx
        File? risk_alleles_genotype_txt = RiskAllelesWorkflowAlias.genotype_txt
        # annotated VCFs
        File qceval_vcf_gz = QCEvalTask.output_vcf_gz
        File? gnomad_vcf_gz = AnnotateGnomadTask.output_vcf_gz
        # GDH export file
        File gdh_export_file = GDHIngestAndFilterTask.matching_variants
        # GDH summary file
        File gdh_summary_file = BgwgsGDHSummaryTask.gdh_summary_file
        File gdh_summary_xlsx = BgwgsGDHSummaryTask.gdh_summary_xlsx
        # IGV report
        File? igv_report = IgvReportFromParsedFASTOutputTask.igv_report
        # GDH Parsed output and NVA report
        File? gdh_parsed_output = GDHOutputParserTask.parsed_report
        File nva_report = GDHOutputParserTask.nva_report
    }
}

task BgwgsGDHSummaryTask {
    input {
        String subject_id
        String sample_id
        String filter_name_or_code
    }

    command <<<
        summary_file="~{subject_id}_~{sample_id}_~{filter_name_or_code}_GDH_summary.txt"

        printf "PM_number\tFilter_Code\n" >> "${summary_file}"
        printf "%s\t%s\n" "~{subject_id}" "~{filter_name_or_code}" >> "${summary_file}"

        pip install openpyxl
        python -c 'import openpyxl
wb = openpyxl.Workbook()
ws = wb.active
ws.cell(row=1, column=1).value = "PM_number"
ws.cell(row=1, column=2).value = "Filter_Code"
ws.cell(row=2, column=1).value = "~{subject_id}"
ws.cell(row=2, column=2).value = "~{filter_name_or_code}"
wb.save("~{subject_id}_~{sample_id}_~{filter_name_or_code}_GDH_summary.xlsx")'
    >>>

    runtime {
        docker: "python:3.10"
    }

    output {
        File gdh_summary_file = subject_id + "_" + sample_id + "_" + filter_name_or_code + "_GDH_summary.txt"
        File gdh_summary_xlsx = subject_id + "_" + sample_id + "_" + filter_name_or_code + "_GDH_summary.xlsx"
    }
}
