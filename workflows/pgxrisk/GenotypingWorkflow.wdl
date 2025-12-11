version 1.0
import "./GATKWorkflow.wdl" as gatk

workflow GenotypingWorkflow {

    input {
        // Inputs for HaplotypeCaller
        File input_cram
        File input_crai
        String sample_id
        String accession_id
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
    String out_prefix = accession_id + "_" + sample_id
    String gvcf = out_path + out_prefix + ".g.vcf"
    String all_calls_vcf = out_path + out_prefix + '.allcalls.vcf'
    String annotated_vcf = out_path + out_prefix + '.annotated.vcf'
    String genotyping_txt = out_path + out_prefix + '_genotyping.txt'

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
            mgbpmbiofx_docker_image = genotyping_docker_image,
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

    call ConvertVcfToTxtTask {
        input:
            input_vcf = AddAnnotationsTask.annotated_vcf_file,
            output_name = genotyping_txt,
            gatk_path = gatk_path,
            genotyping_docker_image = genotyping_docker_image
    }

    output {
        File annotated_vcf_file = AddAnnotationsTask.annotated_vcf_file
        File genotyping_text_report = ConvertVcfToTxtTask.text_report
    }
}

task AddAnnotationsTask {
    input {
        File all_calls_vcf_file
        String annotated_vcf
        String out_path
        String genotyping_docker_image
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
        memory: "24 GB"
    }

    output {
        File annotated_vcf_file = "~{annotated_vcf}"
    }
}

task ConvertVcfToTxtTask {
    input {
        File input_vcf
        String output_name
        String gatk_path
        String genotyping_docker_image
    }

    command <<<
        set -euxo pipefail

        ~{gatk_path} VariantsToTable \
            -V ~{input_vcf} \
            -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
            -O ~{output_name}.txt
    >>>

    runtime {
        docker: "~{genotyping_docker_image}"
        memory: "8 GB"
    }

    output {
        File text_report = "~{output_name}.txt"
    }
}