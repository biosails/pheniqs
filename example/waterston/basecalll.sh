#!/bin/zsh

time bcl2fastq \
--adapter-stringency 1 \
--minimum-trimmed-read-length 1 \
--mask-short-adapter-reads 1 \
-R 170115_NS500488_0313_AHGGKLBGX2 \
-o FASTQ/170115_NS500488_0313_AHGGKLBGX2 \
--sample-sheet FASTQ/170115_NS500488_0313_AHGGKLBGX2_sample_sheet.csv \
--create-fastq-for-index-reads \
--with-failed-reads \
-p 24

time bcl2fastq \
--adapter-stringency 1 \
--minimum-trimmed-read-length 1 \
--mask-short-adapter-reads 1 \
-R 170120_NS500488_0317_AHCFJ2BGX2 \
-o FASTQ/170120_NS500488_0317_AHCFJ2BGX2 \
--sample-sheet FASTQ/170120_NS500488_0317_AHCFJ2BGX2_sample_sheet.csv \
--create-fastq-for-index-reads \
--with-failed-reads \
-p 24

time bcl2fastq \
--adapter-stringency 1 \
--minimum-trimmed-read-length 1 \
--mask-short-adapter-reads 1 \
-R 170123_NS500272_0238_AHGGN5BGX2 \
-o FASTQ/170123_NS500272_0238_AHGGN5BGX2 \
--sample-sheet FASTQ/170123_NS500272_0238_AHGGN5BGX2_sample_sheet.csv \
--create-fastq-for-index-reads \
--with-failed-reads \
-p 24

time bcl2fastq \
--adapter-stringency 1 \
--minimum-trimmed-read-length 1 \
--mask-short-adapter-reads 1 \
-R 170124_NS500272_0239_AH7HYGBGX2 \
-o FASTQ/170124_NS500272_0239_AH7HYGBGX2 \
--sample-sheet FASTQ/170124_NS500272_0239_AH7HYGBGX2_sample_sheet.csv \
--create-fastq-for-index-reads \
--with-failed-reads \
-p 24

time bcl2fastq \
--adapter-stringency 1 \
--minimum-trimmed-read-length 1 \
--mask-short-adapter-reads 1 \
-R 170125_NS500488_0321_AHHJYJBGX2 \
-o FASTQ/170125_NS500488_0321_AHHJYJBGX2 \
--sample-sheet FASTQ/170125_NS500488_0321_AHHJYJBGX2_sample_sheet.csv \
--create-fastq-for-index-reads \
--with-failed-reads \
-p 24
