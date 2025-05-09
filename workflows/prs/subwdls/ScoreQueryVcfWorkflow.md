# RunPRS Workflow

Computes raw and population-adjusted polygenic risk scores (PRS) for
the samples mentioned in a query VCF file.

## Input Parameters

| Type    | Name                      | Req'd | Description | Default Value |
| :---    | :---                      | :---  | :---        | :--- |
| File    | query_vcf                 | Yes   | gz-compressed VCF file of the samples to be scored | |
| File    | adjustment_model_manifest | Yes   | description of the computed model's parameters in JSON format, as that produced by the MakeAdjustmentModel workflow | |
| String  | name                      | Yes   | identifying label for the query dataset | |
| Boolean | norename                  | No    | If `true`, do not run `HelperTasks.RenameChromosomesInVcf` on `query_vcf` | false |

## Output Parameters

| Type | Name            | When   | Description |
| :--- | :---            | :---   | :---        |
| File | raw_scores      | Always | unadjusted PRS for the samples in the query VCF |
| File | adjusted_scores | Always | adjusted PRS for the samples in the query VCF |
| File | pc_projection   | Always | table of reference PC coordinates for the query samples (TSV) |
| File | pc_plot         | Always | plot of the projection of the query samples' PRS onto the space of the first two reference PC (PNG) |
| File | kept_variants   | Always | sites used for training the model |
