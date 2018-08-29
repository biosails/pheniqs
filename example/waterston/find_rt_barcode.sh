#!/usr/bin/env zsh

for fastq in Undetermined_*_R1_001.fastq.gz;
do echo $fastq;
gzcat Undetermined_S0_L001_R1_001.fastq.gz|grep -E '^[ATCGN]{18}$'|sort|uniq -c|sort -rn > ${fastq/.fastq.gz/.uniq}
done
