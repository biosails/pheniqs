#!/usr/bin/env zsh

# lane 4 has 190244417 PF reads

# filter non PF reads and interleave into CRAM
# for lane in {1..8}
# do
# echo "pheniqs mux \
# --base-input /volume/albireo/waverly/CBJLFACXX \
# --base-output /volume/albireo/canal/CBJLFACXX \
# --input Lane${lane}_S${lane}_L00${lane}_I1_001.fastq.gz \
# --input Lane${lane}_S${lane}_L00${lane}_R1_001.fastq.gz \
# --input Lane${lane}_S${lane}_L00${lane}_R2_001.fastq.gz \
# --output CBJLFACXX_l0${lane}.cram \
# --report CBJLFACXX_l0${lane}_report.json \
# 2>&1 > CBJLFACXX_filter.log"
# done

# demultiplex each lane
for lane in {1..8}
do
pheniqs mux --no-output-npf \
--config ~/code/biosails/pheniqs/example/CBJLFACXX/CBJLFACXX_lane${lane}.json \
--output CBJLFACXX_l0${lane}_demultiplexed.bam \
--report CBJLFACXX_l0${lane}_report.json \
2>&1 > CBJLFACXX_demux.log
done
