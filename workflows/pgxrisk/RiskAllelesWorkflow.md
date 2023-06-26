# Risk Alleles Workflow
Runs the workflow to automate discovery, annotation, and report creation for Risk variants for exomes and genomes.

## Input Parameters
| Type | Name | Req'd | Description | 
| :--- | :--- | :---: | :--- |
| File | input_cram | Yes | Sample CRAM/BAM file to run the workflow on |
| File | input_crai | Yes | Corresponding CRAI/BAI file to the CRAM/BAM file |
| String | sample_id | Yes | Full sample identifier |
| String | accession_id | Yes | Sample accession (PM#) |
| String | test_code | Yes | Workflow test-code |
| File | reference_fasta | Yes | HG38 Reference FASTA file |
| File | reference_fasta_fai | Yes | HG38 Reference FASTA index file |
| File | reference_dict | Yes | HG38 Reference dict file |
| File | roi_bed | Yes | Bed file with the workflow region of interest |
| File | dbsnp | Yes | HG38 DBSNP vcf file |
| File | dbsnp_vcf_index | Yes | HG38 DBSNP vcf index file |
| File | workflow_fileset | Yes | Tar file with files pertaining to workflow test-code |
| String | mgbpmbiofx_docker_image | Yes | Name of docker image for the PGx workflow |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | risk_report | Always | Risk Excel report |
| File | genotype_xlsx | Always | Genotype Excel file |
| File | genotype_txt | Always | Genotype text-format file |