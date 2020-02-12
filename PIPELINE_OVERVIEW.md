# Sentieon CCDG & Variant calling

This set of workflows provides a Google Cloud/Sentieon-optimized method of running a CCDG-compliant variant calling pipeline. It is optimized for large sample numbers.
The input for the pipeline is a series of fastq files for individual samples, and the output is a single joint-called VCF per chromosome including all samples. VCFs are split by chromosome to make them less unwieldy with large sample numbers. If you are running this pipeline using smaller sample numbers (<3000), some hard-coded disk sizes may likely be lowered.

Contact heather@dnastack.com for help and clarification.



## Workflows

You will need to set up a Sentieon licence server; the external IP and port for this server should be passed in to workflows 01-05 as the value for `sentieon_licence_server`. (e.g. 111.0.0.1:3000). If you provide this value, you need not provide the other optional Sentieon licence values (`sentieon_auth_mech`, `sentieon_licence_key`).



### [01_fastq_to_gvcf_cram.wdl](https://github.com/DNAstack/mssng-db6-ccdg-alignment-and-calling/blob/master/01_fastq_to_gvcf_cram.wdl)

CCDG-compliant pipeline that converts fastq files into CRAM for long-term storage, and produes a per-sample gVCF file. This file will be combined with other samples in later steps for joint genotyping.
Select `run_genotyper = false` to produce gVCFs for each sample (rather than genotyping single samples).



### [02_joint_genotype_by_region.wdl](https://github.com/DNAstack/mssng-db6-joint-genotype-by-region/blob/1.0.0/02_joint_genotype_by_region.wdl)

Joint genotype samples on a per-region basis. Use `docker/sentieon-bcftools\:201808.06/scripts/generate_shards.sh` to generate region files from a reference genome index or dict using a specified region size (e.g. 50 million base pairs will split the genome into 65 shards). This workflow should then be run once per region. `gvcf_URLs` is a file specifying the gs:// bucket locations of each of the gvcfs output by step 1 (one per line).



### [03_joint_genotype_merge_by_chromosome.wdl](https://github.com/DNAstack/mssng-db6-joint-genotype-merge-by-chr/blob/1.0.0/03_joint_genotype_merge_by_chromosome.wdl)

Merge shards corresponding to a single chromosome, then extract just those regions. `gvcf_URLs` is still required only for sample information. `region` should be a .bed file containing only single chromosome regions (e.g. `chr1.bed` contains `chr1	1	248956422`) (these files for GRCh38 are included in `bed_region_files/`). This will output one `GVCFtyper_main` and one `GVCFtyper_file` file per chromosome. the `GVCFtyper_main` file contains only columns 1-9 of a valid VCF; the `GVCFtyper_file` file contains all sample calls. This allows VQSR to be performed on the (much much smaller) `GVCFtyper_main` file, rather than on the entire final VCF.



### [04_VQSR.wdl (optional)](https://github.com/DNAstack/mssng-db6-vqsr/blob/1.0.0/04_VQSR.wdl)

Perform VQSR on the `GVCFtyper_main` files - they will be combined into a single file, VQSR performed, and then re-split by chromosome. All chromosome `GVCFtyper_main` files should be input here, and one `GVCFtyper_main.recal` recalibrated main file per chromosome will be output. Also input the `chr__.bed` region files from the previous step.



### [05_extract_valid_chromosome_vcfs.wdl](https://github.com/DNAstack/mssng-db6-exract-valid-chromosome-vcfs/blob/1.0.0/05_extract_valid_chromosome_vcfs.wdl)

Produce valid VCFs, one per chromosome. The `vcf_samples_input` is the `GVCFtyper_file` files from step 03 and the `vcf_main` input is either the `GVCFtyper_main` files also from step 03 or the `GVCFtyper_main.recal` recalibrated files from step 04 (depending on whether or not you want VQSR to be performed).
`region` is a string containing the chromosome name only (e.g. chr1). 
`all_samples` is a TSV file with a single line with all sample names separated by tabs - the final VCF will only have these samples in it.
`is_recalibrated` only influences the name of the output file, and refers to whether you're extracting valid VCFs from the base `GVCFtyper_main` file or the recalibrated file.
This will output one final, valid VCF per chromosome.



### [06_variant_metrics.wdl (optional)](https://github.com/DNAstack/mssng-db6-variant-metrics/blob/1.0.0/06_variant_metrics.wdl)

Run once per VCF file output from step 05.
`region` is a string containing chromosome name (same as step05, e.g. `chr10`).
