# Create PRS Variant Weights Files Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | pgs_id | Yes | PGS ID from the PGS Catalog that correspond to a harmonized weights file | |
| File | b38_lookup_file | Yes | GRCh38 lookup file for creating weights files | |
| File | b37_lookup_file | Yes | GRCh37 lookup file for creating weights files | |
| File | chain_file | Yes | Chain file for liftover tool to convert GRCh37 weights files to GRCh38 | |
| String | prs_docker_image | Yes | Docker image for PRS tools | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20260923" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | output_weights | Always | Reformatted GRCh38 weights file |
