version 1.0
import "./GATKWorkflow.wdl" as gatk

workflow GenotypingWorkflow {

    input {
        File input_cram
        File input_crai
        String sample_id
        String accession_id
        String test_code
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        File dbsnp_vcf_index
        String gatk_path = "/gatk/gatk"
        String mgbpmbiofx_docker_image
    }

    String out_path = "outputs/"
    String out_prefix = accession_id + "_" + sample_id + "_" + test_code
    String gvcf = out_path + out_prefix + ".g.vcf"
    String all_calls_vcf = out_path + out_prefix + '.allcalls.vcf'
    String annotated_vcf = out_path + out_prefix + '.annotated.vcf'

    call gatk.HaplotypeCallerTask as HaplotypeCallerTask {
        input:
            input_cram = input_cram,
            input_crai = input_crai,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            roi_bed = roi_bed,
            dbsnp = dbsnp,
            dbsnp_vcf_index = dbsnp_vcf_index,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            gatk_path = gatk_path,
            out_path = out_path,
            gvcf = gvcf,
            all_calls_vcf = all_calls_vcf
    }

    call AddAnnotationsTask {
        input:
            all_calls_vcf_file = HaplotypeCallerTask.all_calls_vcf_file,
            annotated_vcf = annotated_vcf,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
    }

    output {
        File annotated_vcf_file = AddAnnotationsTask.annotated_vcf_file
    }
}

task AddAnnotationsTask {
    input {
        File all_calls_vcf_file
        String annotated_vcf
        String mgbpmbiofx_docker_image
    }

    command <<<
        set -euxo pipefail

        python /biofx-qceval/bin/annotate_with_caller.py \
            ~{all_calls_vcf_file} \
            ~{annotated_vcf}
        
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "24 GB"
    }

    output {
        File annotated_vcf_file = "~{annotated_vcf}"
    }
}