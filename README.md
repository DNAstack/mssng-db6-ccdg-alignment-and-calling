### MSSNG DB6 Sentieon CCDG & Variant calling
CCDG-compliant pipeline that converts fastq files into CRAM for long-term storage, and produes a per-sample gVCF file. This file will be combined with other samples in later steps for joint genotyping.
Select `run_genotyper = false` to produce gVCFs for each sample (rather than genotyping single samples).