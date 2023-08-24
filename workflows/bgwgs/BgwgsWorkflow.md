# BGWGS Workflow
The BGWGS (bigwigs) Workflow starts with a single CRAM or BAM file and provides variant calling,
filtration, annotation, coverage, pharmacogenetic, risk and reporting capabilities.

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230719" |
| String | bcftools_docker_image | No | The name of the bcftools Docker image for VCF annotation | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17" |
| String | subject_id | Yes | The subject id associated with the data | |
| String | sample_id | Yes | The sample id associated with the data | |
| String | sample_data_location | Yes | The cloud storage URL where the sample source data is located | |
| Boolean | fetch_cram | No | Whether or not to fetch the CRAM (primarily for testing) | true |
| Array[String] | fetch_cram_filter_keys | No | The list of strings that must appear in the CRAM file path | [subject_id, sample_id] |
| Array[FileMatcher]? | fetch_cram_file_matchers | No | The list of file matchers to use for CRAM file fetching | |
| Boolean | fetch_bam | No | Whether or not to fetch the BAM (primarily for testing) | true |
| Array[String] | fetch_bam_filter_keys | No | The list of strings that must appear in the BAM file path | [subject_id, sample_id] |
| Array[FileMatcher]? | fetch_bam_file_matchers | No | The list of file matching rules for BAM fetching | |
| Array[String] | fetch_vcf_filter_keys | No | The list of strings that must appear in the VCF file path | [subject_id, sample_id] |
| Array[FileMatcher]? | fetch_vcf_file_matchers | No | The list of file matching rules for VCF fetching | |
| Boolean | fetch_files_verbose | No | If true, generate verbose output from file fetch tasks | false |
| Boolean | do_variant_calling | No | Whether or not to generate a VCF by calling variants from BAM or CRAM; if false, the VCF is fetched from the sample_data_location | true |
| Boolean | do_coverage | No | Whether or not to run depth of coverage analysis | true |
| Boolean | do_pgx | No | Whether or not to generate pharmacogenomics report | true |
| Boolean | do_risk_alleles | No | Whether or not to generate risk alleles report | true |
| Boolean | do_gnomad | No | Whether or not to annotate/load gnomAD data | true |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| File | ref_dict | Yes | The genome reference dict file | |
| File | ref_fasta | Yes | The genome reference fasta file | |
| File | ref_fasta_index | Yes | The genome reference fasta index file | |
| File | dbsnp_vcf | No | dbSNP VCF file that matches the genome reference; required for PGx and Risk Alleles | |
| File | dbsnp_vcf_index | No | dbSNP VCF file index | |
| File | scattered_calling_intervals_list | No | File containing list of scattered calling interval files for Haplotype Caller | |
| File | cov_roi_bed | No | The BED that defines the region of interest for coverage analysis; see DepthOfCoverage.md for details | "gs://lmm-reference-data/roi/clinicalROI_b38_padding2_withCHR_withM_2022-03-11.bed" |
| Array[RoiAndRefGeneFilePair] | cov_roi_genes | No | List of ROI and ref gene file pairs; see DepthOfCoverage.md for details | See Below |
| File | cov_gene_names | No | Tab-delimited file of gene information; see DepthOfCoverage.md for details | "gs://lmm-reference-data/roi/HGNC_genenames_05272022.txt" |
| String | cov_docker_image | No | The name of the coverage Docker image to run the coverage summary task | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/coverage:20230630" |
| String | gatk3_docker_image | No | The name of the GATK3 Docker image for coverage analysis | "broadinstitute/gatk3:3.7-0" |
| String | pgx_test_code | No | Test code that defines which pharmacogenomics report to generate | "lmPGX-pnlC_L" |
| String | pgx_docker_image | No | The name of the Docker image to generate the pharmacogenomics report | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20230711" |
| File | pgx_workflow_fileset | No | Tar file containing the pharmacogenomics reference data to generate the report | "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_files-20220118.tar" |
| File | pgx_roi_bed | No | BED file that defines the genomic regions to include in the pharmacogenomics analysis | "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_genotyping-chr-20220118.bed" |
| String | risk_alleles_test_code | No | Test code that defines which risk alleles report to generate | "lmRISK-pnlB_L" |
| String | risk_alleles_docker_image | No | The name of the Docker image to generate the risk alleles report | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20230724" |
| File | risk_alleles_workflow_fileset | No | Tar file containing the risk alleles reference data to generate the report | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar" |
| File | risk_alleles_roi_bed | No | BED file that defines the genomic regions to include in the risk alleles analysis | "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed" |
| File | target_roi_bed | No | The BED that defines the target region of interest for annotation and filtration | "gs://lmm-reference-data/roi/targetROI_hg38_2023_08_24_withCHR.bed" |
| File | alamut_db | No | The database file for Alamut batch | "gs://lmm-reference-data/annotation/alamut/alamut_db-1.5-2022.01.12.db" |
| File | alamut_fields_tsv | No | The file that defines how the Alamut output is transformed back to a VCF | |
| String | alamut_db_name | No | The database name for the Alamut batch ini file | "alamut_db" |
| String | alamut_server | No | The server name for the Alamut batch ini file | "a-ht-na.interactive-biosoftware.com" |
| String | alamut_port | No | The server port for the Alamut batch ini file | "80" |
| String | alamut_user_secret_name | No | The GCP secret name that contains the user stanza for the Alamut batch ini file | "alamut-batch-ini-user" |
| Int | alamut_queue_limit | No | The maximum number of concurrent Alamut batch processes permitted | 4 |
| String | alamut_queue_folder | No | The shared storage location for Alamut concurrency management | "gs://biofx-task-queue/alamut" |
| String | alamut_docker_image | No | The name of the Alamut Docker image using to run Alamut Batch task | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/alamut:20230630" |
| Boolean | alamut_save_working_files | No | Whether or not to retain intermediate Alamut Batch task files | false |
| String | alamut_anno_src_id | No | When removing already annotated variants prior to Alamut annotation, the annotation source id to query for | "228" |
| String | alamut_anno_min_age | No | When removing already annotated variants prior to Alamut annotation, the annotation minimum timestamp to query for (ISO8601 duration) | "P6M" |
| String | qceval_project_type | No | The type of rules to apply for the QC evaluation task, one of "WGS", "WGS_DRAGEN", "WES" or "NONE" | "WGS" |
| String | qceval_docker_image | No | The name of the Docker image to run the QC evaluation task | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20230630" |
| File | gnomad_coverage_file | No | The gnomad coverage data file | "gs://lmm-reference-data/annotation/gnomad/genomes.r3.0.1.coverage_targetROI-filtered.dedup.txt.gz" |
| File | gnomad_coverage_file_idx | No | The gnomad coverage data index file | "gs://lmm-reference-data/annotation/gnomad/genomes.r3.0.1.coverage_targetROI-filtered.dedup.txt.gz.tbi" |
| Array[String] | gnomad_headers | No | List of VCF headers to add when annotating VCF with gnomad coverage data.  For example, ##INFO=<ID=DP_gnomadG,Number=1,Type=Float,Description="Read depth of GnomAD Genome"> | [ "##INFO=<ID=DP_gnomadG,Number=1,Type=Float,Description=\"Read depth of GnomAD Genome\">" ] |
| String | gnomad_column_list | No | The column list to pass to bcftools annotate for gnomad coverage annotation | "CHROM,POS,INFO/DP_gnomadG" |
| Boolean | has_haploid_sites | No | If true, modify the VCF file headers prior to FAST load to work around lack of support Number=G fields and haploid sites | false |
| String | sample_data_load_config_name | No | The FAST load configuration name for the sample data VCF | "Sample_VCF_PPM_Eval" |
| String | gnomad_data_load_config_name | No | The FAST load configuration name for the gnomad coverage VCF | "Coverage" |
| String | alamut_data_load_config_name | No | The FAST load configuration name for the Alamut annotated VCF | "Alamut" |
| Array[String] | fast_annotated_sample_data_regions | No | The list of regions to include in the FAST annotated sample data; each element is a "name:applyMask" pair |  |
| Array[String] | fast_annotated_sample_data_scripts | No | The list of custom scripts to run on the FAST annotated sample data after creation | |
| String | fast_annotated_sample_data_saved_filter_name | No | The saved filter to apply to the FAST annotated sample data | |
| String | igvreport_docker_image | No | The name of the Docker image to run the IGV report task | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511" |
| String | fast_parser_image | No | The name of the Docker image to run the FAST output parser task | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20230630" |
| File | gil_transcript_exon_count | No | A tab delimited file of transcript id and exon count | "gs://lmm-reference-data/annotation/gil_lmm/transcript_exonNum.txt" |
| String | fast_parser_sample_type | No | The sample type flag for the FAST output parser: S for single-sample Exome or M for multi-sample Exome or B for batch/Biobank or N for NVA-Lite | "S" |

`cov_roi_genes` parameter default:
```json
[
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
```

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | cov_wgs_sample_summary | If coverage analysis is enabled | Coverage metrics summary for all bases |
| File | cov_wgs_sample_statistics | If coverage analysis is enabled | Coverage metrics statistics for all bases |
| File | cov_roi_sample_interval_summary | If coverage analysis is enabled | Coverage metrics interval summary for region of interest |
| File | cov_roi_sample_interval_statistics | If coverage analysis is enabled   Coverage metrics interval statistics for region of interest |
| File | cov_roi_sample_statistics | If coverage analysis is enabled | Coverage metrics statistics for region of interest |
| File | cov_roi_sample_summary | If coverage analysis is enabled | Coverage metrics summary for region of interest |
| File | cov_roi_sample_cumulative_coverage_counts | If coverage analysis is enabled | Coverage cumulative counts for the region of interest |
| File | cov_roi_sample_cumulative_coverage_proportions | If coverage analysis is enabled | Coverage cumulative proportions for the region of interest |
| File | cov_mt_summary | If coverage analysis is enabled | Mitochondrial gene coverage summary |
| File | cov_gene_summary | If coverage analysis is enabled | Gene coverage summary |
| File | cov_gene_summary_unknown | If coverage analysis is enabled | Unknown entries from the gene coverage summary |
| File | cov_gene_summary_entrez | If coverage analysis is enabled | Gene coverage summary enriched with Entrez IDs |
| File | vcf | If variant calling is run | Variants called from BAM/CRAM |
| File | pgx_CPIC_report | If PGx is enabled | CPIC pharmacogenomics report |
| File | pgx_FDA_report | If PGx is enabled | FDA pharmacogenomics report |
| File | pgx_genotype_xlsx | If PGx is enabled | Full list of pharmacogenomics genotypes in XLSX format |
| File | pgx_genotype_txt | If PGx is enabled | Full list of pharmacogenomics genotypes in TSV format |
| File | risk_alleles_report | If risk alleles is enabled | Risk alleles report |
| File | risk_alleles_genotype_xlsx | If risk alleles is enabled | Full list of risk allele genotypes in XLSX format |
| File | risk_alleles_genotype_txt | If risk alleles is enabled | Full list of risk allele genotypes in TSV format |
| File | target_vcf_gz | Always | VCF file filtered to the target region of interest |
| File | alamut_vcf_gz | Always | Target VCF file with Alamut annotations |
| File | qceval_vcf_gz | Always | Target VCF file annotated with QC Evaluation |
| File | gnomad_vcf_gz | If gnomad coverage is enabled | Target VCF annotated with gnomad coverage data |
| File | fast_export_file | Always | Tab-delimited export of annotated sample data from FAST |
| File | fast_summary_file | Always | Summary of FAST processing parameters |
| File | igv_report | Always | HTML-based IGV report file |
| File | fast_parsed_output | Always | Parsed FAST export |
| File | nva_report | Always | NVA report Excel document |