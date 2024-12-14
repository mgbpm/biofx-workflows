version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/HaplotypeCallerGvcfGATK4.wdl"
import "../pgxrisk/PGxWorkflow.wdl"
import "../pgxrisk/RiskAllelesWorkflow.wdl"
import "../../steps/DepthOfCoverage.wdl"
import "../../steps/AlamutBatch.wdl"
import "../../steps/VCFUtils.wdl"
import "../../steps/FASTUtils.wdl"
import "../../steps/QCEval.wdl"
import "../../steps/IgvReport.wdl"
import "../../steps/FASTOutputParser.wdl"
import "../../steps/Utilities.wdl"

workflow BgwgsWorkflow {
    input {
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230828"
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
        # vcf filter inputs
        File target_roi_bed = "gs://lmm-reference-data/roi/targetROI_hg38_2023_08_24_withCHR.bed"
        # alamut inputs
        File alamut_db = "gs://lmm-reference-data/annotation/alamut/alamut_db-1.5-2022.01.12.db"
        File? alamut_fields_tsv
        String alamut_db_name = "alamut_db"
        String alamut_server = "a-ht-na.interactive-biosoftware.com"
        String alamut_port = "80"
        String alamut_user_secret_name = "alamut-batch-ini-user"
        Int alamut_queue_limit = 4
        String alamut_queue_folder = "gs://biofx-task-queue/alamut"
        Int alamut_queue_wait_limit_hrs = 16
        String alamut_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/alamut:20230630"
        Boolean alamut_save_working_files = false
        String alamut_anno_src_id = "228"
        String alamut_anno_min_age = "P6M"
        # qceval inputs
        String qceval_project_type = "WGS"
        String qceval_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20230630"
        # gnomad annotation inputs
        File gnomad_coverage_file = "gs://lmm-reference-data/annotation/gnomad/genomes.r3.0.1.coverage_targetROI-filtered.dedup.txt.gz"
        File gnomad_coverage_file_idx = "gs://lmm-reference-data/annotation/gnomad/genomes.r3.0.1.coverage_targetROI-filtered.dedup.txt.gz.tbi"
        Array[String] gnomad_headers = [ "##INFO=<ID=DP_gnomadG,Number=1,Type=Float,Description=\"Read depth of GnomAD Genome\">" ]
        String gnomad_column_list = "CHROM,POS,INFO/DP_gnomadG"
        # FAST loading inputs
        Boolean has_haploid_sites = false
        String sample_data_load_config_name = "Sample_VCF_PPM_Eval"
        String gnomad_data_load_config_name = "Coverage"
        String alamut_data_load_config_name = "Alamut"
        Array[String] fast_annotated_sample_data_regions
        Array[String]? fast_annotated_sample_data_scripts
        String fast_annotated_sample_data_saved_filter_name
        Int fast_data_load_wait_interval_secs = 300
        Int fast_data_load_wait_max_intervals = 144
        Int fast_adi_wait_interval_secs = 600
        Int fast_adi_wait_max_intervals = 144
        # Reporting steps
        String igvreport_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511"
        String fast_parser_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20230920"
        File gil_transcript_exon_count = "gs://lmm-reference-data/annotation/gil_lmm/transcript_exonNum.txt"
        String fast_parser_sample_type = "S"
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

    # Filter called VCF to target region
    call VCFUtils.FilterVCFWithBEDTask {
        input:
            input_vcf = sample_vcf,
            input_bed = target_roi_bed,
            output_basename = subject_id + "_" + sample_id + ".target",
            docker_image = "~{bcftools_docker_image}"
    }

    # Annotate target VCF with QCEval INFO field
    call QCEval.QCEvalTask {
        input:
            input_vcf = FilterVCFWithBEDTask.output_vcf_gz,
            project_type = qceval_project_type,
            output_basename = subject_id + "_" + sample_id + ".qceval",
            docker_image = qceval_docker_image
    }

    # Pre-filter the VCF file to input to Alamut to remove
    #  variants that already have Alamut annotations
    call FASTUtils.FASTRemoveAlreadyAnnotatedFromVCFTask {
        input:
            input_vcf = FilterVCFWithBEDTask.output_vcf_gz,
            output_basename = subject_id + "_" + sample_id + ".alamutprefilter",
            reference_build = reference_build,
            annotation_source_id = alamut_anno_src_id,
            annotation_min_age = alamut_anno_min_age,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }

    # Annotate target VCF with Alamut
    call AlamutBatch.AlamutBatchTask {
        input:
            input_vcf = FASTRemoveAlreadyAnnotatedFromVCFTask.output_vcf_gz,
            output_basename = subject_id + "_" + sample_id + ".alamut",
            alamut_db = alamut_db,
            alamut_fields_tsv = alamut_fields_tsv,
            alamut_db_name = alamut_db_name,
            alamut_server = alamut_server,
            alamut_port = alamut_port,
            reference_build = reference_build,
            gcp_project_id = gcp_project_id,
            alamut_user_secret_name = alamut_user_secret_name,
            alamut_queue_limit = alamut_queue_limit,
            alamut_queue_folder = alamut_queue_folder,
            alamut_queue_wait_limit_hrs = alamut_queue_wait_limit_hrs,
            docker_image = alamut_docker_image,
            output_working_files = alamut_save_working_files
    }

    # Annotated target VCF with Gnomad coverage
    if (do_gnomad) {
        call VCFUtils.AnnotateVCFTask as AnnotateGnomadTask {
            input:
                input_vcf = FilterVCFWithBEDTask.output_vcf_gz,
                output_basename = subject_id + "_" + sample_id + ".gnomad",
                annotations_file = gnomad_coverage_file,
                annotations_idx_file = gnomad_coverage_file_idx,
                headers_file = write_lines(gnomad_headers),
                column_list = gnomad_column_list,
                docker_image = "~{bcftools_docker_image}"
        }
    }

    # Load sample data to FAST
    call FASTUtils.FASTDataLoadTask as QCEvalLoadTask {
        input:
            reference_build = reference_build,
            vcf_file = QCEvalTask.output_vcf_gz,
            has_haploid_sites = has_haploid_sites,
            sample_data_name = subject_id + "_" + sample_id,
            lab_batch_name = sample_id,
            data_load_config_name = sample_data_load_config_name,
            data_load_target = "SAMPLE_DATA",
            annotation_record_ts = "now",
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }

    # Load Alamut annotations to FAST
    call FASTUtils.FASTDataLoadTask as AlamutLoadTask {
        input:
            reference_build = reference_build,
            vcf_file = AlamutBatchTask.output_vcf_gz,
            has_haploid_sites = has_haploid_sites,
            data_load_config_name = alamut_data_load_config_name,
            data_load_target = "ANNOTATION_DATA",
            merge_strategy = "MERGE",
            annotation_record_ts = "now",
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }

    # Load Gnomad annotations to FAST
    if (defined(AnnotateGnomadTask.output_vcf_gz)) {
        call FASTUtils.FASTDataLoadTask as GnomadLoadTask {
            input:
                reference_build = reference_build,
                vcf_file = select_first([AnnotateGnomadTask.output_vcf_gz]),
                has_haploid_sites = has_haploid_sites,
                data_load_config_name = gnomad_data_load_config_name,
                data_load_target = "ANNOTATION_DATA",
                merge_strategy = "MERGE",
                annotation_record_ts = "now",
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }
    }

    # Wait for up to 12 hours data loads to complete
    call FASTUtils.FASTWaitForDataLoadsTask {
        input:
            data_load_ids = select_all([QCEvalLoadTask.data_load_id, AlamutLoadTask.data_load_id, GnomadLoadTask.data_load_id]),
            check_interval_secs = fast_data_load_wait_interval_secs,
            max_checks = fast_data_load_wait_max_intervals,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    if (FASTWaitForDataLoadsTask.wait_result.total != FASTWaitForDataLoadsTask.wait_result.completed) {
        call Utilities.FailTask as DataLoadTimeout {
            input:
                error_message = "One or more FAST data loads did not complete within allotted time"
        }
    }
    if (FASTWaitForDataLoadsTask.wait_result.total != FASTWaitForDataLoadsTask.wait_result.succeeded) {
        call Utilities.FailTask as DataLoadFailure {
            input:
                error_message = "One or more FAST data loads did not succeed"
        }
    }

    # Create annotated sample data
    if (FASTWaitForDataLoadsTask.wait_result.total == FASTWaitForDataLoadsTask.wait_result.succeeded) {

        # Wait up to 24 hours for annotation data initialization to complete
        call FASTUtils.FASTWaitForADITask {
            input:
                check_interval_secs = fast_adi_wait_interval_secs,
                max_checks = fast_adi_wait_max_intervals,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }
        if (FASTWaitForADITask.num_pending >= 1) {
            call Utilities.FailTask as ADINotComplete {
                input:
                    error_message = "FAST annotation data initialization did not complete within allotted time"
            }
        }
        if (FASTWaitForADITask.num_pending <= 0) {
            # Create annotated sample data
            call FASTUtils.FASTCreateAnnotatedSampleDataTask {
                input:
                    annotated_sample_data_name = subject_id + "_" + sample_id + "_" + fast_annotated_sample_data_saved_filter_name,
                    sample_data_names_and_labels = [subject_id + "_" + sample_id],
                    region_names_and_masks = fast_annotated_sample_data_regions,
                    scripts = fast_annotated_sample_data_scripts,
                    saved_filter_name = fast_annotated_sample_data_saved_filter_name,
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name,
                    docker_image = orchutils_docker_image
            }

            # Wait up to 6 hours for annotated sample data to complete
            call FASTUtils.FASTWaitForAnnotatedSampleDataTask {
                input:
                    annotated_sample_data_name = FASTCreateAnnotatedSampleDataTask.annotated_sample_data_name_output,
                    check_interval_secs = 300,
                    max_checks = 72,
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name,
                    docker_image = orchutils_docker_image
            }
            if (FASTWaitForAnnotatedSampleDataTask.wait_result.total != FASTWaitForAnnotatedSampleDataTask.wait_result.completed) {
                call Utilities.FailTask as AnnotatedSampleTimeout {
                    input:
                        error_message = "FAST annotated sample data creation did not complete within allotted time"
                }
            }
            if (FASTWaitForAnnotatedSampleDataTask.wait_result.total != FASTWaitForAnnotatedSampleDataTask.wait_result.succeeded) {
                call Utilities.FailTask as AnnotatedSampleFailure {
                    input:
                        error_message = "FAST annotated sample data creation did not succeed"
                }
            }

            if (FASTWaitForAnnotatedSampleDataTask.wait_result.total == FASTWaitForAnnotatedSampleDataTask.wait_result.succeeded) {
                call FASTUtils.FASTExportAnnotatedSampleDataTask {
                    input:
                        annotated_sample_data_name = FASTCreateAnnotatedSampleDataTask.annotated_sample_data_name_output,
                        format = "TXT",
                        output_basename = FASTCreateAnnotatedSampleDataTask.annotated_sample_data_name_output + ".fastexport",
                        gcp_project_id = gcp_project_id,
                        workspace_name = workspace_name,
                        docker_image = orchutils_docker_image
                }
                call FASTOutputParser.FASTOutputParserTask {
                    input:
                        fast_output_file = FASTExportAnnotatedSampleDataTask.output_file,
                        sample_type = fast_parser_sample_type,
                        reference_build = reference_build,
                        oms_query = "Y",
                        transcript_exonNum = gil_transcript_exon_count,
                        report_basename = FASTCreateAnnotatedSampleDataTask.annotated_sample_data_name_output,
                        gcp_project_id = gcp_project_id,
                        workspace_name = workspace_name,
                        fast_parser_image = fast_parser_image
                }
                if (defined(FASTOutputParserTask.parsed_report)) {
                    call IgvReport.IgvReportFromParsedFASTOutputTask {
                        input:
                            bam_cram = sample_bam,
                            bai_crai = sample_bai,
                            parsed_fast_output = select_first([FASTOutputParserTask.parsed_report]),
                            output_basename = FASTCreateAnnotatedSampleDataTask.annotated_sample_data_name_output + ".igvreport",
                            ref_fasta = ref_fasta,
                            ref_fasta_index = ref_fasta_index,
                            track_files = igv_track_files,
                            track_index_files = igv_track_index_files,
                            docker_image = igvreport_docker_image
                    }
                }
                call BgwgsFASTSummaryTask {
                    input:
                        subject_id = subject_id,
                        sample_id = sample_id,
                        fast_annotated_sample_data_regions = fast_annotated_sample_data_regions,
                        fast_annotated_sample_data_saved_filter_name = fast_annotated_sample_data_saved_filter_name
                }
            }
        }
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
        # filtered VCF
        File target_vcf_gz = FilterVCFWithBEDTask.output_vcf_gz
        # annotated VCFs
        File alamut_vcf_gz = AlamutBatchTask.output_vcf_gz
        File qceval_vcf_gz = QCEvalTask.output_vcf_gz
        File? gnomad_vcf_gz = AnnotateGnomadTask.output_vcf_gz
        # FAST export file
        File? fast_export_file = FASTExportAnnotatedSampleDataTask.output_file
        # FAST summary file
        File? fast_summary_file = BgwgsFASTSummaryTask.fast_summary_file
        File? fast_summary_xlsx = BgwgsFASTSummaryTask.fast_summary_xlsx
        # IGV report
        File? igv_report = IgvReportFromParsedFASTOutputTask.igv_report
        # FAST Parsed output and NVA report
        File? fast_parsed_output = FASTOutputParserTask.parsed_report
        File? nva_report = FASTOutputParserTask.nva_report
    }
}

task BgwgsFASTSummaryTask {
    input {
        String subject_id
        String sample_id
        Array[String]? fast_annotated_sample_data_regions
        String fast_annotated_sample_data_saved_filter_name
    }

    command <<<
        summary_file="~{subject_id}_~{sample_id}_~{fast_annotated_sample_data_saved_filter_name}_FAST_summary.txt"

        printf "PM_number\tSaved_filter\tRegions\n" >> "${summary_file}"
        printf "%s\t%s\t%s\n" "~{subject_id}" "~{fast_annotated_sample_data_saved_filter_name}" "~{sep=',' fast_annotated_sample_data_regions}" >> "${summary_file}"

        pip install openpyxl
        python -c 'import openpyxl
wb = openpyxl.Workbook()
ws = wb.active
ws.cell(row=1, column=1).value = "PM_number"
ws.cell(row=1, column=2).value = "Saved_filter"
ws.cell(row=1, column=3).value = "Regions"
ws.cell(row=2, column=1).value = "~{subject_id}"
ws.cell(row=2, column=2).value = "~{fast_annotated_sample_data_saved_filter_name}"
ws.cell(row=2, column=3).value = "~{sep=',' fast_annotated_sample_data_regions}"
wb.save("~{subject_id}_~{sample_id}_~{fast_annotated_sample_data_saved_filter_name}_FAST_summary.xlsx")'
    >>>

    runtime {
        docker: "python:3.10"
    }

    output {
        File fast_summary_file = subject_id + "_" + sample_id + "_" + fast_annotated_sample_data_saved_filter_name + "_FAST_summary.txt"
        File fast_summary_xlsx = subject_id + "_" + sample_id + "_" + fast_annotated_sample_data_saved_filter_name + "_FAST_summary.xlsx"
    }
}