## VCF Call Biobank Crams

A workflow to process the Biobank CRAMs stored in Wasabi into haplotype called vcf.

Input .crams have a single identifier.
Output .vcf will have a second identifier added

The CopyFilesTask will be modified to rename the .cram on copy to Terra

eg: 10003797.cram -> 10003797_f639b25627af47bd9042a173f363d49b.cram

All Samples will be loaded into Sample Tracker and will include 2 identifiers.

