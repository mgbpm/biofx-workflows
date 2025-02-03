# SubsetShards Workflow
Workflow to subset the shard files under a specified input directory tree, and (optionally) concatenate the resulting subsets.

### Input Parameters
* String sourcebase - required - the root of the directory tree where the shard files to subset are located; these files are identified by two conditions: (1) having a `.vcf.gz` extension; and (2) having the component `/shards/` in their paths
* String? shardsbase - optional - the root of the directory tree where FIXME
* String targetbase - required - the root of the directory tree where the requested outputs will be stored
* File shardsmapfile - required - JSON file describing the shards to subset and the desired outputs (see full description below)
* Boolean noconcat - optional - if true, the generated subset shards will not be concatenated; defaults to false
* Boolean? nodelete - optional - if true, the generated shards will not be deleted at the end of the workflow; FIXME
* File? regions - optional - file specifying regions to extract from the input shards; see below for default behavior
* File? samples - optional - file specifying samples to extract from the input shards; see below for default behavior
* File? script - optional - script to carry out the subsetting operation on each shard; if omitted, the VM's `subsetvcf.sh` script will be used; (FIXME: reference to this script?)
* Int nbatches - optional - number of batches to split the subsetting operation in; defaults to 500
* String docker_image - optional - docker image to use for the workflow; defaults to "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
* Int ConcatenateShards.memory - optional - gigabytes of memory for ConcatenateShards task; defaults to 16

### Output Parameters
None.  (The workflow's results are saved under `targetbase`.)

## Description

### The shardsmap file

The shards map file is a JSON file with the following structure

    [
        {
            "shards": [
                "<relpath-to-shard-0-0>",
                "<relpath-to-shard-0-1>",
                ...
                "<relpath-to-shard-0-M0>"
            ],
            "vcf": "<relpath-to-result-0>"
        },
        {
            "shards": [
                ...
            ],
            "vcf": "<relpath-to-result-1>"
        },

        ...
        {
            "shards": [
                ...
            ],
            "vcf": "<relpath-to-result-N>"
        }
    ]

In words, the JSON specifies a an array of dictionaries, each of which has two members `shards` and `vcf`.  The value of each `shards` member is an array of shard paths relative to the root specified by the workflow input `sourcebase`.  The value of each `vcf` member is a path relative to the root specified by the workflow input `targetbase`.  It specifies the URL for the concatenation of the subset shards generated from the shards specified in the corresponding `shards` member.  Following the example above, the shards at

    <sourcebase>/<relpath-to-shard-0-0>
    <sourcebase>/<relpath-to-shard-0-1>
    ...
    <sourcebase>/<relpath-to-shard-0-M0>

...will be subset, and the resulting subset shards will be concatenated, in the ordered specified in the `shards` member, and the result saved as

    <targetbase>/<relpath-to-result-N>
