# GenotypingWorkflow

This workflow performs variant genotyping for exome or genome samples, producing a fully annotated VCF file.

## Input Parameters

| Type   | Name                   | Req'd | Description |
| :----- | :--------------------- | :---: | :---------- |
| File   | input_cram             | Yes   | Sample CRAM/BAM file to run the workflow on |
| File   | input_crai             | Yes   | Corresponding CRAI/BAI file to the CRAM/BAM file |
| String | sample_id              | Yes   | Full sample identifier |
| String | accession_id           | Yes   | Sample accession (PM#) |
| String | test_code              | Yes   | Workflow test-code |
| File   | reference_fasta        | Yes   | HG38 Reference FASTA file |
| File   | reference_fasta_fai    | Yes   | HG38 Reference FASTA index file |
| File   | reference_dict         | Yes   | HG38 Reference dict file |
| File   | roi_bed                | Yes   | Bed file with the workflow region of interest |
| File   | dbsnp                  | Yes   | HG38 DBSNP VCF file |
| File   | dbsnp_vcf_index        | Yes   | HG38 DBSNP VCF index file |
| String | gatk_path              | No    | Path to GATK executable (default: /gatk/gatk) |
| String | mgbpmbiofx_docker_image| Yes   | Name of docker image for the workflow |

## Output Parameters

| Type | Name                | When   | Description |
| :--- | :------------------ | :----- | :---------- |
| File | annotated_vcf_file  | Always | VCF file with genotypes and caller annotation |

## Workflow Steps

1. **HaplotypeCallerTask**  
   Runs GATK HaplotypeCaller to generate a GVCF and an all-calls VCF.

2. **AddAnnotationsTask**  
   Annotates the VCF with caller/version information using `annotate_with_caller.py`.

## Example Output File Names

- `<accession_id>_<sample_id>_<test_code>.g.vcf`
- `<accession_id>_<sample_id>_<test_code>.allcalls.vcf`
- `<accession_id>_<sample_id>_<test_code>.allcalls.vcf.annotated.vcf`

---