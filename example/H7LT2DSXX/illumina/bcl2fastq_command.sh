#!/bin/zsh

bcl2fastq \
--runfolder-dir 181014_A00534_0024_AH7LT2DSXX \
--sample-sheet no_barcode_samplesheet.csv \
--create-fastq-for-index-reads \
--fastq-compression-level 4 \
--output-dir H7LT2DSXX \
--adapter-stringency 0 \
--minimum-trimmed-read-length 0 \
--mask-short-adapter-reads 0
# --no-bgzf-compression
