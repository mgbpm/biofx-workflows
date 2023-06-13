version 1.0

workflow AlamutBatchWorkflow {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".alamut"
        File alamut_db
        File? alamut_fields_tsv
        String alamut_db_name = "alamut_db"
        String alamut_server = "a-ht-na.interactive-biosoftware.com"
        String alamut_port = "80"
        String reference_build = "GRCh38"
        String gcp_project_id
        String alamut_user_secret_name = "alamut-batch-ini-user"
        Int alamut_queue_limit = 4
        String alamut_queue_folder = "gs://biofx-task-queue/alamut"
        String docker_image
        Int disk_size = 150
        Boolean output_working_files = false
    }

    call AlamutBatchTask {
        input:
            input_vcf = input_vcf,
            output_basename = output_basename,
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
            docker_image = docker_image,
            disk_size = disk_size,
            output_working_files = output_working_files
    }

    output {
        File output_vcf_gz = AlamutBatchTask.output_vcf_gz
    }
}

task AlamutBatchTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".alamut"
        File alamut_db
        File? alamut_fields_tsv
        String alamut_db_name = "alamut_db"
        String alamut_server = "a-ht-na.interactive-biosoftware.com"
        String alamut_port = "80"
        String reference_build = "GRCh38"
        String gcp_project_id
        String alamut_user_secret_name = "alamut-batch-ini-user"
        Int alamut_queue_limit = 4
        String alamut_queue_folder = "gs://biofx-task-queue/alamut"
        String docker_image
        Int disk_size = 150
        Boolean output_working_files = false
    }

    command <<<
        set -euxo pipefail

        # make directory for working files
        mkdir workdir

        # ensure we are starting with an uncompressed VCF
        bcftools view --output-type v "~{input_vcf}" > workdir/input.vcf

        # prepare alamut ini and db
        # [Database]
        # File=/data/pcpgm/share/annot_db/alamut/alamut_db-1.5-2022.01.12.db
        # Name=alamut_db
        # [Network]
        # IBS\Server=a-ht-na.interactive-biosoftware.com
        # IBS\Port=80
        # [User]
        # Institution=PCPGM_HT
        # LicenceKey=53421643
        # LicenseAccepted=true
        # User=PCPGM
        ALAMUT_BATCH_EXE=$(which alamut-batch)
        ALAMUT_BATCH_DIR=$(dirname "$ALAMUT_BATCH_EXE")
        ALAMUT_BATCH_INI="${ALAMUT_BATCH_DIR}/alamut-batch.ini"
        echo "[Database]" > "${ALAMUT_BATCH_INI}"
        echo "File=~{alamut_db}" >> "${ALAMUT_BATCH_INI}"
        echo "Name=~{alamut_db_name}" >> "${ALAMUT_BATCH_INI}"
        echo "" >> "${ALAMUT_BATCH_INI}"
        echo "[Network]" >> "${ALAMUT_BATCH_INI}"
        echo "IBS\Server=~{alamut_server}" >> "${ALAMUT_BATCH_INI}"
        echo "IBS\Port=~{alamut_port}" >> "${ALAMUT_BATCH_INI}"
        echo "" >> "${ALAMUT_BATCH_INI}"

        gcloud --project="~{gcp_project_id}" secrets versions access "latest" --secret=~{alamut_user_secret_name} --out-file="${ALAMUT_BATCH_INI}.user"
        cat "${ALAMUT_BATCH_INI}.user" >> "${ALAMUT_BATCH_INI}"
        echo "" >> "${ALAMUT_BATCH_INI}"

        # prepare input VCF for alamut
        cat workdir/input.vcf | sed "s/AD,Number=./AD,Number=R/" | \
            $MGBPMBIOFXPATH/biofx-alamut/bin/vcf_from_gvcf.pl 0 | \
            $MGBPMBIOFXPATH/biofx-alamut/bin/vcf_wide_to_tall.pl | \
            $MGBPMBIOFXPATH/biofx-alamut/bin/vcf_ala_id.pl > workdir/alamut.input.vcf

        # concurrency management, wait until we can get slot in the queue
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w global -r "~{alamut_queue_folder}"
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/enter_task_queue.py --queue-folder "~{alamut_queue_folder}" \
            --entry-details "source = AlamutBatch.wdl, input vcf = ~{input_vcf}" \
            --queue-limit ~{alamut_queue_limit} \
            --wait --wait-ttl-hrs 8 --wait-interval-mins 5 | tee queue_entry_id
        QUEUE_ENTRY_ID=$(cat queue_entry_id)

        # run alamut batch process (but don't abort script on failure to cleanup queue)
        set +e
        alamut-batch --in workdir/alamut.input.vcf --ann workdir/alamut.output.ann.tsv --unann workdir/alamut.output.unann.tsv \
            --assbly ~{reference_build} --processes 1 --alltrans --ssIntronicRange 2 --nomispred --outputannonly
        ALA_RETVAL=$?

        # remove the queue entry
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/remove_task_queue_entry.py --queue-folder "~{alamut_queue_folder}" --entry-id ${QUEUE_ENTRY_ID}

        # restore the fast fail option
        set -e

        # alamut error handling
        [ ${ALA_RETVAL} -ne 0 ] && exit ${ALA_RETVAL}

        # post process alamut output from TSV to VCF
        #  using the field mapping to determine which fields are common and which are transcript specific
        sed -i "/LRG/d" workdir/alamut.output.ann.tsv
        ALA_FIELDS_TSV="~{alamut_fields_tsv}"
        if [ -z "${ALA_FIELDS_TSV}" ]
        then
            ALA_FIELDS_TSV="$MGBPMBIOFXPATH/biofx-alamut/config/alamut_fields_v1.2.0.tsv"
        fi
        bcftools view -h workdir/alamut.input.vcf | grep -i ^##contig > workdir/vcf.contig.headers.txt
        $MGBPMBIOFXPATH/biofx-alamut/bin/parse_alamut_fieldbyname_v2.pl workdir/alamut.output.ann.tsv "${ALA_FIELDS_TSV}" workdir/vcf.contig.headers.txt > workdir/alamut.output.ann.vcf
        bcftools sort --output-type z workdir/alamut.output.ann.vcf > "~{output_basename}.vcf.gz"
    >>>

    runtime {
        docker: "~{docker_image}"
        memory: "8GB"
        disks: "local-disk " + disk_size + " SSD"
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
        Array[File] working_files = if output_working_files then glob("workdir/*") else []
    }
}