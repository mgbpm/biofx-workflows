# Alamut Batch Step
The Alamut batch step runs the Alamut batch program to annotate a VCF file with
reference data from the provided Alamut database.

## Input Parameters
* File input_vcf - required - the VCF file to annotate
* String output_basename - optional - the output file basename, defaults to the input_vcf name with the file extension(s) dropped
* File alamut_db - required - the Alamut database file
* File alamut_fields_tsv - optional - file that defines the rules for transforming the Alamut TSV output into a VCF file,
  defaults to biofx-alamut/config/alamut_fields_v1.2.0.tsv if not specified
* String alamut_db_name - optional - the name of the Alamut database for the configuration file, defaults to "alamut_db"
* String alamut_server - optional - the Alamut license server hostname, defaults to "a-ht-na.interactive-biosoftware.com"
* String alamut_port - optional - the Alamut license server port, defaults to "80"
* String reference_build - optional - the reference build for the VCF, defaults to "GRCh38"
* String gcp_project_id - required - the GCP project that stores the Alamut user configuration secret
* String alamut_user_secret_name - optional - the name of the secret that contains the [User] stanza for the configuration,
  defaults to "alamut-batch-ini-user"
* Int alamut_queue_limit - optional - the maximum number of concurrent Alamut tasks to allow, defaults to 4
* String alamut_queue_folder - optional - the Alamut concurrency queue folder (cloud storage URI), 
  defaults to "gs://biofx-task-queue/alamut"
* Int alamut_queue_wait_limit_hrs - optional - the maximum number of hours to wait for a queue slot before failing
* String alamut_docker_image - required - the name of the Alamut docker image to run

## Output Parameters
* File output_vcf_gz - the annotated VCF file

## Docker Image Requirements
* Executables
  * alamut-batch
  * bcftools - for VCF file manipulation
  * python3 - Python 3.x for task concurrency management
  * perl - for Alamut pre- and post-processing scripts
  * gcloud - for Alamut configuration secrets
* MGBPM Biofx Code
  * biofx-alamut - helper scripts for Alamut pre- and post-processing
  * biofx-orchestration-utils - task concurrency management scripts