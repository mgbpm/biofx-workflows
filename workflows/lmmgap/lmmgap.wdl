version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/Utilities.wdl"

workflow LmmgapWorkflow {
    input {
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        String orchutils_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
        String lmmgap_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/lmmgap:20231206"
        String batch_id
        String project_file_location
        Array[String] subject_ids
        Array[String] sample_ids
        Array[String] sample_genders
        Array[String] sample_data_locations
        File genome_reference_fa
        File genome_reference_fai
        File genome_reference_dict
        String genome_reference_name = "hg19"
        File reference_sample_project_file = "gs://lmm-reference-data/lmmgap/27202778065-8_ReferenceProjectDetailReport.csv"
        String reference_sample_gender = "2"
        String reference_subject_id = "NA12878-A10"
        String reference_sample_id = "204339030100-R01C01"
        File reference_sample_gtc_file = "gs://lmm-reference-data/lmmgap/NA12878-A10_204339030100-R01C01.gtc"
        Float call_rate_cut_off = 0.99
        File bpm_manifest_file = "gs://lmm-reference-data/lmmgap/GDA-8v1-0_A5.bpm"
        File csv_manifest_file = "gs://lmm-reference-data/lmmgap/GDA-8v1-0_A5.csv"
        File excluded_site_file = "gs://lmm-reference-data/lmmgap/GDA-8v1-0_A5.probes-to-exclude.tsv"
        Boolean save_working_files = true
    }

    # Fetch the project details file
    call FileUtils.FetchFilesTask as FetchProjectFile {
        input:
            data_location = sub(project_file_location, "/[^/]+$", ""),
            recursive = true,
            file_match_keys = [ basename(project_file_location) ],
            docker_image = orchutils_docker_image,
            disk_size = 15,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name
    }
    if (length(FetchProjectFile.all_files) <= 0) {
        call Utilities.FailTask as ProjectFileNotFound {
            input:
                error_message = "Project file " + project_file_location + " not found"
        }
    }

    # Fetch the GTC files
    scatter (i in range(length(sample_ids))) {
        call FileUtils.FetchFilesTask as FetchGTCFile {
            input:
                data_location = sample_data_locations[i],
                recursive = true,
                file_types = [ "gtc" ],
                file_match_keys = [ subject_ids[i], sample_ids[i] ],
                docker_image = orchutils_docker_image,
                disk_size = 15,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name
        }
        if (length(FetchGTCFile.all_files) <= 0) {
            call Utilities.FailTask as GTCFileNotFound {
                input:
                    error_message = "GTC file for sample " + sample_ids[i] + " not found in " + sample_data_locations[i]
            }
        }
    }

    # Map sample tracker sex codes to plink
    call Utilities.MapSampleTrackerSexToPlinkSexTask {
        input:
            sample_tracker_sexes = sample_genders
    }

    # Run lmm gap
    call LmmgapTask {
        input:
            run_id = batch_id,
            project_details_file = FetchProjectFile.all_files[0],
            subject_ids = subject_ids,
            sample_ids = sample_ids, 
            sample_genders = MapSampleTrackerSexToPlinkSexTask.plink_sexes,
            gtc_files = flatten(FetchGTCFile.all_files),
            reference_sample_project_file = reference_sample_project_file,
            reference_subject_id = reference_subject_id,
            reference_sample_id = reference_sample_id,
            reference_sample_gender = reference_sample_gender,
            reference_sample_gtc_file = reference_sample_gtc_file,
            call_rate_cut_off = call_rate_cut_off,
            genome_reference_fa = genome_reference_fa,
            genome_reference_fai = genome_reference_fai,
            genome_reference_dict = genome_reference_dict,
            genome_reference_name = genome_reference_name,
            bpm_manifest_file_path = bpm_manifest_file,
            csv_manifest_file_path = csv_manifest_file,
            excluded_site_file_path = excluded_site_file,
            docker_image = lmmgap_docker_image,
            save_working_files = save_working_files
    }

    output {
        File final_report = LmmgapTask.final_report
        File master_file = LmmgapTask.master_file
        Array[File] sample_vcfs = LmmgapTask.sample_vcfs
        Array[File] sample_vcf_idxs  = LmmgapTask.sample_vcf_idxs
        File merged_vcf = LmmgapTask.merged_vcf
        File merged_bed = LmmgapTask.merged_bed
        File merged_bed_bim = LmmgapTask.merged_bed_bim
        File merged_bed_map = LmmgapTask.merged_bed_map
        File merged_bed_hh =  LmmgapTask.merged_bed_hh
        File merged_bed_ped = LmmgapTask.merged_bed_ped
        File merged_bed_fam = LmmgapTask.merged_bed_fam
        File split_bed = LmmgapTask.split_bed
        File split_bed_hh = LmmgapTask.split_bed_hh
        File split_bed_bim = LmmgapTask.split_bed_bim
        File split_bed_fam = LmmgapTask.split_bed_fam
        File genome_qc = LmmgapTask.genome_qc
        File genome_qc_hh = LmmgapTask.genome_qc_hh
        File sexcheck_qc = LmmgapTask.sexcheck_qc
        File sexcheck_qc_hh = LmmgapTask.sexcheck_qc_hh
        File missing_qc_imiss = LmmgapTask.missing_qc_imiss
        File missing_qc_lmiss = LmmgapTask.missing_qc_lmiss
        File missing_qc_hh = LmmgapTask.missing_qc_hh
        Array[File] working_gtc_to_vcf_files = LmmgapTask.working_gtc_to_vcf_files
        Array[File] working_fixhdr_vcf_files = LmmgapTask.working_fixhdr_vcf_files
        Array[File] working_sorted_vcf_files = LmmgapTask.working_sorted_vcf_files
        Array[File] working_fixhet_vcf_files = LmmgapTask.working_fixhet_vcf_files
        Array[File] working_log_files = LmmgapTask.working_log_files
    }
}   
task LmmgapTask {
    input {
        String run_id
        File project_details_file
        Array[String] subject_ids
        Array[String] sample_ids
        Array[String] sample_genders
        Array[File] gtc_files
        File reference_sample_project_file
        String reference_subject_id
        String reference_sample_id
        String reference_sample_gender
        File reference_sample_gtc_file
        Float call_rate_cut_off
        File genome_reference_fa
        File genome_reference_fai
        File genome_reference_dict
        String genome_reference_name
        File bpm_manifest_file_path
        File csv_manifest_file_path
        File excluded_site_file_path
        String docker_image
        Int preemptible = 1
        Boolean save_working_files = true
    }

    command <<<
        set -euxo pipefail
        
        # array tools are very fussy about the genome reference files
        #  the extensions must be .fa, .fai and .dict
        #  and symlinks are not allowed
        mkdir genome_reference
        mv "~{genome_reference_fa}" genome_reference/reference.fa
        mv "~{genome_reference_fai}" genome_reference/reference.fa.fai
        mv "~{genome_reference_dict}" genome_reference/reference.dict

        # build the gender file
        subject_ids=("~{sep='" "' subject_ids}")
        sample_ids=("~{sep='" "' sample_ids}")
        genders=("~{sep='" "' sample_genders}")
        len=${#subject_ids[@]}
        for (( i=0; i<len; i++ ))
        do
            gender_code=$(echo "${genders[$i]}" | cut -c1)
            echo -e "${subject_ids[$i]}\t${sample_ids[$i]}\t${gender_code}" >> genders.tsv
        done
        echo -e "~{reference_subject_id}\t~{reference_sample_id}\t~{reference_sample_gender}" >> genders.tsv

        mkdir output
        ${MGBPMBIOFXPATH}/biofx-lmmgap/bin/pipeline.py \
            --output-directory output \
            --gtc-file-paths "~{sep=',' gtc_files}" \
            --excluded-site-list-path "~{excluded_site_file_path}" \
            --csv-manifest-path "~{csv_manifest_file_path}" \
            --bpm-manifest-path "~{bpm_manifest_file_path}" \
            --genome-reference-path genome_reference/reference.fa \
            --gender-id-file "genders.tsv" \
            --project-details-file "~{project_details_file}" \
            --call-rate-cut-off "~{call_rate_cut_off}" \
            --reference-genome-name "~{genome_reference_name}" \
            --reference-sample-project-file "~{reference_sample_project_file}" \
            --reference-sample-path "~{reference_sample_gtc_file}" \
            --run-id "~{run_id}"
    >>>

    runtime {
        docker: "~{docker_image}"
        # for each GTC file about 1-1.5 GB of disk space is required
        disks: "local-disk " + (ceil(length(gtc_files) * 1.5) + 10) + " SSD"
        memory: "16GB"
        preemptible: preemptible
    }

    output {
        File final_report = "output/report/~{run_id}.final_report.tsv"
        File master_file = "output/project/~{run_id}_master.tsv"
        Array[File] sample_vcfs = glob("output/collapse_vcf/*.vcf.gz")
        Array[File] sample_vcf_idxs  = glob("output/collapse_vcf/*.vcf.gz.tbi")
        File merged_vcf =       "output/merged_vcf/~{run_id}_merged.vcf.gz"
        File merged_bed =       "output/merged_bed/~{run_id}_merged.bed"
        File merged_bed_bim =   "output/merged_bed/~{run_id}_merged.bim"
        File merged_bed_map =   "output/merged_bed/~{run_id}_merged.map"
        File merged_bed_hh =    "output/merged_bed/~{run_id}_merged.hh"
        File merged_bed_ped =   "output/merged_bed/~{run_id}_merged.ped"
        File merged_bed_fam =   "output/merged_bed/~{run_id}_merged.fam"
        File split_bed =        "output/split_bed/~{run_id}_split.bed"
        File split_bed_hh =     "output/split_bed/~{run_id}_split.hh"
        File split_bed_bim =    "output/split_bed/~{run_id}_split.bim"
        File split_bed_fam =    "output/split_bed/~{run_id}_split.fam"
        File genome_qc =        "output/qc/genome/~{run_id}_genome.genome"
        File genome_qc_hh =     "output/qc/genome/~{run_id}_genome.hh"
        File sexcheck_qc =      "output/qc/sex_check/~{run_id}_sex_check.sexcheck"
        File sexcheck_qc_hh =   "output/qc/sex_check/~{run_id}_sex_check.hh"
        File missing_qc_imiss = "output/qc/missing/~{run_id}_missing.imiss"
        File missing_qc_lmiss = "output/qc/missing/~{run_id}_missing.lmiss"
        File missing_qc_hh =    "output/qc/missing/~{run_id}_missing.hh"
        Array[File] working_gtc_to_vcf_files = if save_working_files then glob("output/gtc_to_vcf/*.vcf") else []
        Array[File] working_fixhdr_vcf_files = if save_working_files then glob("output/fixhdr_vcf/*.vcf") else []
        Array[File] working_sorted_vcf_files = if save_working_files then glob("output/sorted_vcf/*.vcf") else []
        Array[File] working_fixhet_vcf_files = if save_working_files then glob("output/fixhet_vcf/*.vcf") else []
        Array[File] working_log_files = if save_working_files then
            [ 
                "output/split_bed/~{run_id}_split.log", 
                "output/merged_bed/~{run_id}_merged.log", 
                "output/qc/missing/~{run_id}_missing.log", 
                "output/qc/genome/~{run_id}_genome.log", 
                "output/qc/sex_check/~{run_id}_sex_check.log"
            ] 
            else []
    }
}
