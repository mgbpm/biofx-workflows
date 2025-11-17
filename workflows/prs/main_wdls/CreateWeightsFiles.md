# Create PRS Variant Weights Files Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[String] | pgs_ids | Yes | PGS IDs from the PGS Catalog that correspond to weights files | |
| File | b38_lookup | Yes | GRCh38 lookup file for creating weights files | |
| File | b37_lookup | Yes | GRCh37 lookup file for creating weights files | |
| File | chain_file | Yes | Chain file for liftover tool to convert GRCh37 weights files to GRCh38 | |
| String | prs_docker_image | Yes | Docker image for PRS tools | |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
