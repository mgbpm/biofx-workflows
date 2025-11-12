# PRS Structs for Tasks

## SelfExclusiveSites

| Type | Name | Req'd | Description |
| :--- | :--- | :---: | :--- |
| File | sites | Yes | Sites must have three columns: ID, chrom, pos |
| Int | maxAllowed | Yes | Number of sites allowed |

## WeightSet

| Type | Name | Req'd | Description |
| :--- | :--- | :---: | :--- |
| File | linear_weights | Yes | A standard PRS weights file |
| File | interaction_weights | No | An interaction weights file; must have these columns in any order: id_1, id_2, chrom_1, chrom_2, pos_1, pos_2, allele_1, allele_2, weight |
| SelfExclusiveSites| interaction_self_exclusive_sites | No | The interaction term will only be added in no more than selfExclusiveSites.maxAllowed of the effect alleles listed in SelfExclusizeSites.sites is observed |

## ScoringInputs

| Type | Name | Req'd | Description |
| :--- | :--- | :---: | :--- |
| File | variant_weights | Yes | A variant weights file |
| File | training_variants | Yes | Contains the variants used to train the ancestry adjustment model |

## AdjustmentModelData

| Type | Name | Req'd | Description |
| :--- | :--- | :---: | :--- |
| String | condition_code | Yes | Condition/disease related to the model |
| File | parameters | Yes | Fitted params from training the model |
| File | principal_components | Yes | PCs from training the model |
| File | loadings | Yes | PC loadings from training the model |
| File | meansd | Yes | PC meansd file from training the model |
| Array[File] | variant_weights | Yes | Variant weights files used to train the model |
| File | score_weights | No | Weights of variant weights files used to train a PRSmix model; not required for single weights PRS models, but is required for PRSmix models |
| Array[ScoringInputs] | scoring_inputs | Yes | ScoringInputs struct with the variant weights files used to train the model and the training variants output from the model |
| File | pca_variants | Yes | PCA variants used to train the model |
| File | original_pca_variants | Yes | All PCA variants used before subsetting to variants common among the input VCFs, weights files, and PCA file |
| File | query_file | Yes | Query VCF used to train the model |
| Int | base_memory | Yes | Estimated base memory required to process the reference |
