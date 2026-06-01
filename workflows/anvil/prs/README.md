# PRS Mix and Single Scoring

A polygenic risk score (PRS) summarizes the chance of developing a certain health condition based on estimated effects of many genetic variants on an individual’s phenotype. In a standard run of the PRS pipeline, a single variant weights file is used to calculate a polygenic risk score. However, to more accurately determine an individual’s score, a mix score can be calculated using multiple variant weights files and a score weights file.

The PRS pipeline in this repo can calculate mix or single scores for multiple health conditions for many individuals. It will calculate the PRS for each provided variant weights file per health condition of interest. If given multiple variant weights files, it will compute the mix score for each health condition by finding the sum of weighted polygenic risk scores associated with the condition. Then, if desired, the PRS mix or single score can then be adjusted using a previously computed ancestry adjustment model and PCA.

Inputs for the PRS pipeline have strict stipulations in order for the pipeline to run successfully. To faciliatate successful runs of the pipeline, we have created a WDL for input preparation and a WDL for training ancestry adjustment models that are also within this repo.

## Folder Structure

There are three main WDLs in the scripts folder that can be used to run PRS mix and single scoring. Each are described in more detail in their corresponding section of the README. The WDLs are:

1. PrsInputPrep.wdl
2. PrsModel.wdl
3. PrsScoring.wdl

The res folder within this repo contains flowcharts that outline the steps taken in each of the three main PRS WDLs. The flowcharts can also be found in the corresponding WDL's section below.

Example JSON input files for each of the three main PRS WDLs can be found in the examples folder of this repo.

## PrsInputPrep WDL

To run the PRS pipeline, input files must be prepared and cleaned to follow format requirements and other pipeline stipulations. The PRS Input Preparation workflow will scrub input files to only include the intersection of variants found in the variant weights files, PCA variants file, reference VCF, and query VCFs. In addition, it will ensure that naming conventions in the input files match those expected by the pipeline.

Below is a flowchart that outlines the workflow step-by-step:
![PRS Input Prep Flowchart](/res/prsprepinputs_flowchart.png)

Input files for the PRS Input Preparation WDL include:

1. **Weights files** -- Variant weights files for the PRS pipeline must be tab-delimited and contain three columns in the following order: variant ID, effect allele, and effect weight. The variant ID should be comprised of the chromosome, position, allele, and effect allele (ex: 1:927392:A:G). The [PGS Catalog](https://www.pgscatalog.org/) is a good resource for variant weights files, however, the files must be amended to the described format. _For the pipeline to calculate a mix score, all weights files for calculating the mix score must be provided and the PGS ID associated with each weights file must be in the corresponding file name._
2. **PCA variants file** -- The PCA variants file is a file with one column containing varaint IDs in chromosome, position, ref, alt notation.
3. **Query VCF files** -- All VCFs that will be scored for a health condition of interest and therefore with the same condition's ancestry adjustment model. The VCFs can be single-sample or contain multiple samples/individuals.

### PRS Prep Input Parameters

| Type        | Name             | Required | Description | Default Value |
| :---        | :---             | :---     | :---        | :---          |
| Array[File] | weights_files    | Yes      | Array of 3-column TSV files describing variants and their weights; provide multiple weights files if mix scoring is desired | |
| File        | pca_variants     | Yes      | Text file listing the variants for principal component analysis (PCA), one variant per line | |
| String      | workspace        | Yes      | Name of the Terra workspace where the workflow will be run | |
| String      | source           | Yes      | URL to location of reference VCF shards | |
| String      | target           | Yes      | URL to use as temporary work space and to save the generated reference VCF | |
| Int         | nbatches         | No       | Number of parallel jobs for subsetting the reference shards | 500 |
| Boolean     | resuming         | No       | Whether this run is the resumption of an earlier run | false |
| Boolean     | norename         | No       | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| Array[File] | query_vcfs       | Yes      | Array of gz-compressed VCF files of the samples to be scored | |
| String      | prs_docker_image | No       | Docker image equipped with PRS scripts | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515" |

### PRS Prep Output Parameters

| Type        | Name                  | When   | Description |
| :---        | :---                  | :---   | :---        |
| File        | regions               | Always | Regions file used to generate reference VCF |
| File        | kept_pca_variants     | Always | File listing the retained and renamed PCA variants |
| Array[File] | renamed_weights_files | Always | Weights files with renamed variant IDs |
| Array[File] | renamed_query_vcfs    | Always | Query VCFs with renamed variant IDs |
| File        | reference_vcf         | Always | Generated refernce VCF |
| File        | reference_tbi         | Always | Index for generated reference VCF |

## PrsModel WDL

An ancestry adjustment model will be used to adjust a mix or single PRS. Multiple models -- one for each health condition of interest -- can be supplied as input to the PRS workflow as a model manifest, and a PRS will be calculated for each model. The model manifests for the PRS pipeline are JSON files that contain the health condition name, paths to variant weights files, a path to a PC loadings file, paths to training variants files, a path to a PCA variants file, a path to the query VCF used to create the model, a path to a meansd file, a path to parameters for the model, and a path to the score weights file. Models and their associated model manifests can be generated using the Prs Model WDL.

Below is a flowchart that outlines the workflow step-by-step:
![PRS Model Flowchart](/res/prsmodel_flowchart.png)

Inputs for the PRS Model WDL include:

1. **Variant weights files** -- The cleaned and prepped variant weights file from the PRS Input Prep WDL should be used as input for the PRS Model WDL. These weights files will be in the same format as the originals, but will have corrected naming conventions.
2. **PCA variants file** -- The output PCA variants file from the PRS Input Prep WDL should be used as input for the PRS Model WDL. This file will be in the same format as the original, but it will be filtered to the intersection of variants in the variant weights files, query files, and original PCA variants file.
3. **Reference VCF** -- We recommend using a reference such as 1KG for the reference VCF. Note that the reference VCF should not be the same as any query VCF. The reference VCF used should be the output reference VCF from the PRS Input Preparation WDL.
4. **Query file** -- The query file for making an ancestry adjustment model should be one of the output query VCFs from the PRS Input Prep WDL.
5. **Score weights** -- A score weights file must be tab-delimited and contain two columns in the following order: PGS ID and weight. The PGS IDs found in input variant weights filenames will be used to associate a weight to the variant weights file when calculating mix scores. To calculate a single PRS and not a mix score, a file with a single PGS ID and weight of 1 can be used as input.

### PRS Model Input Parameters

| Type        | Name                   | Required | Description | Default Value |
| :---------- | :--------------------- | :------: | :---------- | :------------ |
| File        | condition_name         | Yes      | Code for condition/disease | |
| Array[File] | var_weights            | Yes      | Array of different PGS variant weight files, with 3 columns: variant ID, effect allele, and score; include multiple files for mix scoring and only one file for single scoring | |
| File        | pca_variants           | Yes      | Variants used in PCA projection | |
| File        | reference_vcf          | Yes      | Reference VCF of population for creating/training model | |
| File        | query_file             | Yes      | VCF or TSV to score | |
| File        | score_weights          | No       | Score weights for each PGS ID; required for mix models | |
| Boolean     | norename               | No       | If `true`, do not rename chromosomes to have chr prefix | false |

### PRS Model Output Parameters

| Type        | Name                      | When   | Description |
| :---------- | :------------------------ | :----- | :---------- |
| File        | model_manifest            | Always | JSON of model file paths for either single score model or mix model |

## PrsScoring WDL

Once the inputs have been prepared and any desired ancestry adjustment models have been created, the PRS Scoring WDL can be run. The PRS Scoring workflow will generate a PRS per variant weights file for each of the input model manifests. It will then calculate a mix or single score for each health condition associated, and adjust the resulting score using the ancestry adjustment model and PCA (if desired).

Below is a flowchart that outlines the workflow step-by-step:
![PRS Scoring Flowchart](/res/prsscoring_flowchart.png)

### PRS Scoring Input Parameters

| Type        | Name                | Required | Description | Default Value |
| :---------- | :------------------ | :------: | :---------- | :--- |
| File        | input_vcf           | Yes      | VCF to score | |
| Array[File] | model_manifests     | Yes      | Adjustment model manifest file from MakeAdjustmentModelWorkflow WDL | |
| Boolean     | norename            | No       | If `true`, do not rename chromosomes to have chr prefix | false |
| Boolean     | perform_adjustment  | No       | If `true`, use ancestry adjustment model to adjust scores | true |
| File        | renaming_lookup     | No       | Mapping file for renaming chromosomes | "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv" |
| String      | ubuntu_docker_image | No       | Ubuntu Docker image | "ubuntu:latest" |

### PRS Scoring Output Parameters

| Type        | Name               | When   | Description |
| :---------- | :----------------- | :----- | :---------- |
| Array[File] | prs_raw_scores     | Always | PRS raw scores |
| File        | prs_mix_raw_score  | If a score weights file exists in the model manifest | PRS mix raw scores |
| File        | prs_adjusted_score | If perform_adjustment is `true` | PRS scores or mix scores adjusted with population models |

## References

Original PRS Github: https://github.com/broadinstitute/palantir-workflows/tree/main/PRS
