version 1.0

import "./FileUtils.wdl"

workflow GenotypingBGWGSWorkflow {
    input {
        String sample_id
        String accession_id

        # Inputs for the FetchFilesTask
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String? data_location
        Array[String] fetch_cram_filter_keys = [accession_id, sample_id]
        Boolean fetch_cram = false
        Array[FileMatcher]? fetch_cram_file_matchers
        Boolean fetch_files_verbose = false
        String workspace_name
        Int fetch_disk_size = 75
        String orchutils_docker_image = "us-docker.pkg.dev/mgbpmbiofx/orchestration-utils/orchestration-utils:20250203"

        # Inputs for HaplotypeCaller
        File? input_cram
        File? input_crai
        String test_code
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        File dbsnp_vcf_index
        String gatk_path = "/gatk/gatk"
        String genotyping_docker_image
    }

    String out_path = "outputs/"
    String out_prefix = accession_id + "_" + sample_id + "_" + test_code
    String gvcf = out_path + out_prefix + ".g.vcf"
    String all_calls_vcf = out_path + out_prefix + '.allcalls.vcf'
    String annotated_vcf = out_path + out_prefix + '.annotated.vcf'
    
    if (fetch_cram) {
        call FileUtils.FetchFilesTask as FetchCram {
            input:
                data_location = select_first([data_location, input_cram]),
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

    call HaplotypeCallerTask {
        input:
            input_cram = select_first([FetchCram.cram, input_cram]),
            input_crai = select_first([FetchCram.crai, input_crai]),
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            roi_bed = roi_bed,
            dbsnp = dbsnp,
            dbsnp_vcf_index = dbsnp_vcf_index,
            genotyping_docker_image = genotyping_docker_image,
            gatk_path = gatk_path,
            out_path = out_path,
            gvcf = gvcf,
            all_calls_vcf = all_calls_vcf
    }

    call AddAnnotationsTask {
        input:
            all_calls_vcf_file = HaplotypeCallerTask.all_calls_vcf_file,
            annotated_vcf = annotated_vcf,
            out_path = out_path,
            genotyping_docker_image = genotyping_docker_image
    }

    output {
        File annotated_vcf_file = AddAnnotationsTask.annotated_vcf_file
    }
}

task HaplotypeCallerTask {
    input {
        File input_cram
        File input_crai
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        File dbsnp_vcf_index
        String genotyping_docker_image
        String gatk_path
        String out_path
        String gvcf
        String all_calls_vcf
        Int mem_gb = 24
        Int disk_size = ceil(size(input_cram, "GB") + size(dbsnp, "GB")) + 10
    }
    
    command <<<
        set -euxo pipefail

        mkdir -p ~{out_path}

        ~{gatk_path} --java-options "-Xmx20G" \
        HaplotypeCaller \
        --input ~{input_cram} \
        --output ~{gvcf} \
        --reference ~{reference_fasta} \
        --intervals ~{roi_bed} \
        --dbsnp ~{dbsnp} \
        -mbq 10 \
        -stand-call-conf 30 \
        -A QualByDepth \
        -A FisherStrand \
        -A RMSMappingQuality \
        -A MappingQualityZero \
        -A StrandOddsRatio \
        -A DepthPerAlleleBySample \
        --read-filter MappingQualityReadFilter \
        --minimum-mapping-quality '17' \
        --read-filter MappingQualityNotZeroReadFilter \
        -ERC BP_RESOLUTION
        
        ~{gatk_path} --java-options "-Xmx20G" \
        GenotypeGVCFs \
        --variant ~{gvcf} \
        --output ~{all_calls_vcf} \
        --reference ~{reference_fasta} \
        --intervals ~{roi_bed} \
        --dbsnp ~{dbsnp} \
        -A QualByDepth \
        -A FisherStrand \
        -A RMSMappingQuality \
        -A MappingQualityZero \
        -A StrandOddsRatio \
        -A DepthPerAlleleBySample \
        --read-filter MappingQualityReadFilter \
        --minimum-mapping-quality '17' \
        --read-filter MappingQualityNotZeroReadFilter
        
    >>>

    runtime {
        docker: "~{genotyping_docker_image}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
    }

    output {
        File all_calls_vcf_file = "~{all_calls_vcf}"
        File all_calls_vcf_idx_file = "~{all_calls_vcf}.idx"
        File gvcf_file = "~{gvcf}"
        File gvcf_idx_file = "~{gvcf}.idx"
    }
}

task AddAnnotationsTask {
    input {
        File all_calls_vcf_file
        String annotated_vcf
        String out_path
        String genotyping_docker_image
        Int mem_gb = 24
        Int disk_size = ceil(size(all_calls_vcf_file, "GB")) + 10
    }

    command <<<
        set -euxo pipefail

        mkdir -p ~{out_path}

        python $MGBPMBIOFXPATH/biofx-qceval/bin/annotate_with_caller.py \
            ~{all_calls_vcf_file} \
            ~{annotated_vcf}
        
    >>>

    runtime {
        docker: "~{genotyping_docker_image}"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb} GB"
    }

    output {
        File annotated_vcf_file = "~{annotated_vcf}"
    }
}
