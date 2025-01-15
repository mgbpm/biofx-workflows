version 1.0

import "../../steps/Utilities.wdl"
import "../../steps/FileUtils.wdl"
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/PGx_3.0.1/workflows/pgxrisk/PGxWorkflow.wdl" as PGx_v3
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/Risk_3.0.1/workflows/pgxrisk/RiskAllelesWorkflow.wdl" as Risk_v3

workflow BahrainPGxRiskWorkflow {
    input {
        # Sample data inputs
        Array[String] sample_ids
        Array[String] collaborator_sample_ids
        Array[String] data_bucket
        # Orchutils docker iamge
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20240625"
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Reference genome files
        String reference_build = "GRCh38"
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File dbsnp_vcf
        File dbsnp_vcf_index
        # PGx inputs
        String pgx_test_code = "lmPGX-pnlC_L"
        String pgx_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:dev-v3"
        File pgx_workflow_fileset = "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_files-20220118.tar"
        File pgx_roi_bed = "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_genotyping-chr-20220118.bed"
        # Risk alleles inputs
        String risk_alleles_test_code = "lmRISK-pnlB_L"
        String risk_alleles_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20240129"
        File risk_alleles_workflow_fileset = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar"
        File risk_alleles_roi_bed = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed"
    }

    # Fetch CRAM, VCF, and index files
    scatter (i in range(length(sample_ids))) {
        call FileUtils.FetchFilesTask as FetchScreeningFiles {
            input:
                data_location = data_bucket[i],
                recursive = true,
                file_types = [ "cram" ],
                file_match_keys = [ sample_ids[i] ],
                docker_image = orchutils_docker_image,
                disk_size = 50,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name
        }
        if (!defined(FetchScreeningFiles.cram)) {
            call Utilities.FailTask as CRAMNotFound {
                input:
                    error_message = "CRAM for sample " + sample_ids[i] + " not found in " + data_bucket[i]
            }
        }
        if (!defined(FetchScreeningFiles.crai)) {
            call Utilities.FailTask as CRAMIndexNotFound {
                input:
                    error_message = "CRAM index for sample " + sample_ids[i] + " not found in " + data_bucket[i]
            }
        }
    }
    # Coerce object types to Array[File] for future tasks
    Array[File] fetched_crams = select_all(select_first([FetchScreeningFiles.cram]))
    Array[File] fetched_crai = select_all(select_first([FetchScreeningFiles.crai]))


    # Run PGx & Risk workflows with fetched CRAMs
    scatter (i in range(length(fetched_crams))) {
        call PGx_v3.PGxWorkflow as RunPGx {
            input:
                input_cram = fetched_crams[i],
                input_crai = fetched_crai[i],
                sample_id = sample_ids[i],
                accession_id = collaborator_sample_ids[i],
                test_code = pgx_test_code,
                reference_fasta = ref_fasta,
                reference_fasta_fai = ref_fasta_index,
                reference_dict = select_first([ref_dict]),
                roi_bed = pgx_roi_bed,
                dbsnp = select_first([dbsnp_vcf]),
                dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
                workflow_fileset = pgx_workflow_fileset,
                mgbpmbiofx_docker_image = pgx_docker_image
        }
        call Risk_v3.RiskAllelesWorkflow as RunRisk {
            input:
                input_cram = fetched_crams[i],
                input_crai = fetched_crai[i],
                sample_id = sample_ids[i],
                accession_id = collaborator_sample_ids[i],
                test_code = risk_alleles_test_code,
                reference_fasta = ref_fasta,
                reference_fasta_fai = ref_fasta_index,
                reference_dict = select_first([ref_dict]),
                roi_bed = risk_alleles_roi_bed,
                dbsnp = select_first([dbsnp_vcf]),
                dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
                workflow_fileset = risk_alleles_workflow_fileset,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                mgbpmbiofx_docker_image = risk_alleles_docker_image
        }
    }

    
    output {
        # PGx outputs
        Array[File]? pgx_CPIC_report = RunPGx.CPIC_report
        Array[File]? pgx_FDA_report = RunPGx.FDA_report
        Array[File]? pgx_genotype_xlsx = RunPGx.genotype_xlsx
        Array[File]? pgx_genotype_txt = RunPGx.genotype_txt
        # Risk outputs
        Array[File]? risk_alleles_report = RunRisk.risk_report
        Array[File]? risk_alleles_genotype_xlsx = RunRisk.genotype_xlsx
        Array[File]? risk_alleles_genotype_txt = RunRisk.genotype_txt
    }
}