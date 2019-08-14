#!/usr/bin/env zsh

chmod +x H7LT2DSXX_basecall.sh

time ./H7LT2DSXX_basecall.sh
# ./H7LT2DSXX_basecall.sh  680997.56s user 3131.80s system 2361% cpu 8:02:45.45 total

time pheniqs mux --config H7LT2DSXX_l01_estimate.json
time pheniqs mux --config H7LT2DSXX_l02_estimate.json
time pheniqs mux --config H7LT2DSXX_l03_estimate.json
time pheniqs mux --config H7LT2DSXX_l04_estimate.json
# pheniqs mux --config H7LT2DSXX_l01_estimate.json  73651.91s user 25660.06s system 1567% cpu 1:45:37.11 total
# pheniqs mux --config H7LT2DSXX_l02_estimate.json  74040.02s user 25411.91s system 1560% cpu 1:46:11.54 total
# pheniqs mux --config H7LT2DSXX_l03_estimate.json  75065.86s user 25380.67s system 1584% cpu 1:45:40.76 total
# pheniqs mux --config H7LT2DSXX_l04_estimate.json  73312.76s user 25702.11s system 1608% cpu 1:42:36.82 total

time pheniqs mux --config H7LT2DSXX_l01_adjusted.json
time pheniqs mux --config H7LT2DSXX_l02_adjusted.json
time pheniqs mux --config H7LT2DSXX_l03_adjusted.json
time pheniqs mux --config H7LT2DSXX_l04_adjusted.json

# pheniqs mux --config H7LT2DSXX_l01_adjusted.json  204310.87s user 33165.95s system 1146% cpu 5:45:14.91 total
# pheniqs mux --config H7LT2DSXX_l02_adjusted.json  205686.39s user 33449.11s system 1146% cpu 5:47:36.80 total
# pheniqs mux --config H7LT2DSXX_l03_adjusted.json  202190.24s user 33469.03s system 1184% cpu 5:31:29.01 total
# pheniqs mux --config H7LT2DSXX_l04_adjusted.json  200925.97s user 33061.22s system 1183% cpu 5:29:35.39 total

# total reads = 11578868372
# SS = Split Segment
# IS = Interleaved Segment
# SL = Split Library
# IL = Interleaved Library
#
# tool        input format    input layout    output format   output layout   time        Max Mem
# pheniqs     fastq           SS/IL           null            null            01:25:55    19MB
# pheniqs     fastq           SS/IL           bam             IS/IL           03:10:07    51MB
# pheniqs     fastq           SS/IL           fastq           SS/IL           03:20:46    68MB
# pheniqs     fastq           SS/IL           fastq           SS/SL           03:30:25    990MB
# pheniqs     fastq           SS/IL           cram            IS/IL           04:07:47    402MB
# deml        fastq           SS/IL           fastq           SS/SL           17:47:12    212MB
# basecall    bcl                             fastq           SS/IL           00:46:56    5869MB
# bcl2fastq   bcl                             fastq           SS/SL           01:02:40    22005MB
#
#
#
# index tool        input         output      time        memory
# 1     pheniqs     fastq SS/IL   null        01:25:55    19MB
# 2     pheniqs     fastq SS/IL   bam IS/IL   03:10:07    51MB
# 3     pheniqs     fastq SS/IL   fastq SS/IL 03:20:46    68MB
# 4     pheniqs     fastq SS/IL   fastq SS/SL 03:30:25    990MB
# 5     pheniqs     fastq SS/IL   cram IS/IL  04:07:47    402MB
# 6     deml        fastq SS/IL   fastq SS/SL 17:47:12    212MB
# 7     bcl2fastq   bcl           fastq SS/IL 00:46:56    5869MB
# 8     bcl2fastq   bcl           fastq SS/SL 01:02:40    22005MB


# benchmark on prince

# /home/lg1883/.bin/pheniqs mux --config /scratch/lg1883/H7LT2DSXX/H7LT2DSXX_l01_estimate.json --base-input /scratch/lg1883/H7LT2DSXX/one --base-output /scratch/lg1883/H7LT2DSXX/estimate
#     User time (seconds): 48489.46
#     System time (seconds): 84033.73
#     Percent of CPU this job got: 2570%
#     Elapsed (wall clock) time (h:mm:ss or m:ss): 1:25:55
#     Average shared text size (kbytes): 0
#     Average unshared data size (kbytes): 0
#     Average stack size (kbytes): 0
#     Average total size (kbytes): 0
#     Maximum resident set size (kbytes): 19444
#     Average resident set size (kbytes): 0
#     Major (requiring I/O) page faults: 0
#     Minor (reclaiming a frame) page faults: 5403282
#     Voluntary context switches: 916574893
#     Involuntary context switches: 1440279
#     Swaps: 0
#     File system inputs: 86505040
#     File system outputs: 224
#     Socket messages sent: 0
#     Socket messages received: 0
#     Signals delivered: 0
#     Page size (bytes): 4096
#     Exit status: 0

# /home/lg1883/.bin/pheniqs mux --config /scratch/lg1883/H7LT2DSXX/H7LT2DSXX_l01_bam.json --base-input /scratch/lg1883/H7LT2DSXX/one --base-output /scratch/lg1883/H7LT2DSXX/bam
#     User time (seconds): 192213.17
#     System time (seconds): 59427.54
#     Percent of CPU this job got: 2205%
#     Elapsed (wall clock) time (h:mm:ss or m:ss): 3:10:07
#     Average shared text size (kbytes): 0
#     Average unshared data size (kbytes): 0
#     Average stack size (kbytes): 0
#     Average total size (kbytes): 0
#     Maximum resident set size (kbytes): 52164
#     Average resident set size (kbytes): 0
#     Major (requiring I/O) page faults: 4
#     Minor (reclaiming a frame) page faults: 1061313645
#     Voluntary context switches: 3290052858
#     Involuntary context switches: 57440210
#     Swaps: 0
#     File system inputs: 1008654504
#     File system outputs: 991280920
#     Socket messages sent: 0
#     Socket messages received: 0
#     Signals delivered: 0
#     Page size (bytes): 4096
#     Exit status: 0

# /home/lg1883/.bin/pheniqs mux --config /scratch/lg1883/H7LT2DSXX/H7LT2DSXX_l01_fastq.json --base-input /scratch/lg1883/H7LT2DSXX/one --base-output /scratch/lg1883/H7LT2DSXX/fastq
# 	User time (seconds): 190675.50
# 	System time (seconds): 80039.28
# 	Percent of CPU this job got: 2247%
# 	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:20:46
# 	Average shared text size (kbytes): 0
# 	Average unshared data size (kbytes): 0
# 	Average stack size (kbytes): 0
# 	Average total size (kbytes): 0
# 	Maximum resident set size (kbytes): 69656
# 	Average resident set size (kbytes): 0
# 	Major (requiring I/O) page faults: 4
# 	Minor (reclaiming a frame) page faults: 941148482
# 	Voluntary context switches: 3492233233
# 	Involuntary context switches: 75034191
# 	Swaps: 0
# 	File system inputs: 1008654504
# 	File system outputs: 758065568
# 	Socket messages sent: 0
# 	Socket messages received: 0
# 	Signals delivered: 0
# 	Page size (bytes): 4096
# 	Exit status: 0

# Command being timed: "/home/lg1883/.bin/pheniqs mux --config /scratch/lg1883/H7LT2DSXX/H7LT2DSXX_l01_fastq_split.json --base-input /scratch/lg1883/H7LT2DSXX/one --base-output /scratch/lg1883/H7LT2DSXX/fastq_split --buffer 256"
#     User time (seconds): 245041.06
#     System time (seconds): 41281.90
#     Percent of CPU this job got: 2267%
#     Elapsed (wall clock) time (h:mm:ss or m:ss): 3:30:25
#     Average shared text size (kbytes): 0
#     Average unshared data size (kbytes): 0
#     Average stack size (kbytes): 0
#     Average total size (kbytes): 0
#     Maximum resident set size (kbytes): 1013360
#     Average resident set size (kbytes): 0
#     Major (requiring I/O) page faults: 1
#     Minor (reclaiming a frame) page faults: 1065774279
#     Voluntary context switches: 2226487477
#     Involuntary context switches: 144435930
#     Swaps: 0
#     File system inputs: 1008647360
#     File system outputs: 911891208
#     Socket messages sent: 0
#     Socket messages received: 0
#     Signals delivered: 0
#     Page size (bytes): 4096
#     Exit status: 0

# /home/lg1883/.bin/pheniqs mux --config /scratch/lg1883/H7LT2DSXX/H7LT2DSXX_l01_cram.json --base-input /scratch/lg1883/H7LT2DSXX/one --base-output /scratch/lg1883/H7LT2DSXX/cram
# 	User time (seconds): 147206.18
# 	System time (seconds): 114556.62
# 	Percent of CPU this job got: 1760%
# 	Elapsed (wall clock) time (h:mm:ss or m:ss): 4:07:47
# 	Average shared text size (kbytes): 0
# 	Average unshared data size (kbytes): 0
# 	Average stack size (kbytes): 0
# 	Average total size (kbytes): 0
# 	Maximum resident set size (kbytes): 411760
# 	Average resident set size (kbytes): 0
# 	Major (requiring I/O) page faults: 8
# 	Minor (reclaiming a frame) page faults: 448771056
# 	Voluntary context switches: 3426584394
# 	Involuntary context switches: 8183243
# 	Swaps: 0
# 	File system inputs: 1008649960
# 	File system outputs: 554319464
# 	Socket messages sent: 0
# 	Socket messages received: 0
# 	Signals delivered: 0
# 	Page size (bytes): 4096
# 	Exit status: 0

# /home/lg1883/.bin/deML --index /home/lg1883/H7LT2DSXX/H7LT2DSXX_deml_l01_index.txt --outfile /home/lg1883/H7LT2DSXX/deml/H7LT2DSXX_l01 -f /home/lg1883/H7LT2DSXX/one/H7LT2DSXX_S1_L001_R1_001.fastq.gz -r /home/lg1883/H7LT2DSXX/one/H7LT2DSXX_S1_L001_R2_001.fastq.gz -if1 /home/lg1883/H7LT2DSXX/one/H7LT2DSXX_S1_L001_I1_001.fastq.gz -if2 /home/lg1883/H7LT2DSXX/one/H7LT2DSXX_S1_L001_I2_001.fastq.gz
#     User time (seconds): 63553.39
#     System time (seconds): 440.39
#     Percent of CPU this job got: 99%
#     Elapsed (wall clock) time (h:mm:ss or m:ss): 17:47:12
#     Average shared text size (kbytes): 0
#     Average unshared data size (kbytes): 0
#     Average stack size (kbytes): 0
#     Average total size (kbytes): 0
#     Maximum resident set size (kbytes): 217808
#     Average resident set size (kbytes): 0
#     Major (requiring I/O) page faults: 3
#     Minor (reclaiming a frame) page faults: 5360211
#     Voluntary context switches: 12343
#     Involuntary context switches: 52211
#     Swaps: 0
#     File system inputs: 1030202808
#     File system outputs: 99669832
#     Socket messages sent: 0
#     Socket messages received: 0
#     Signals delivered: 0
#     Page size (bytes): 4096
#     Exit status: 0

# bcl2fastq -l DEBUG --tiles s_1_ --runfolder-dir /scratch/lg1883/H7LT2DSXX/181014_A00534_0024_AH7LT2DSXX --sample-sheet /scratch/lg1883/H7LT2DSXX/H7LT2DSXX_basecall_sample_sheet.csv --create-fastq-for-index-reads --adapter-stringency 0 --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --output-dir /scratch/lg1883/H7LT2DSXX --fastq-compression-level 3
#     User time (seconds): 76895.92
#     System time (seconds): 816.72
#     Percent of CPU this job got: 2759%
#     Elapsed (wall clock) time (h:mm:ss or m:ss): 46:55.78
#     Average shared text size (kbytes): 0
#     Average unshared data size (kbytes): 0
#     Average stack size (kbytes): 0
#     Average total size (kbytes): 0
#     Maximum resident set size (kbytes): 6009452
#     Average resident set size (kbytes): 0
#     Major (requiring I/O) page faults: 14
#     Minor (reclaiming a frame) page faults: 102756941
#     Voluntary context switches: 122688
#     Involuntary context switches: 4065166
#     Swaps: 0
#     File system inputs: 667680328
#     File system outputs: 1007714496
#     Socket messages sent: 0
#     Socket messages received: 0
#     Signals delivered: 0
#     Page size (bytes): 4096
#     Exit status: 0

# bcl2fastq -l DEBUG --tiles s_1_ --runfolder-dir /scratch/lg1883/H7LT2DSXX/181014_A00534_0024_AH7LT2DSXX --sample-sheet /scratch/lg1883/H7LT2DSXX/181014_A00534_0024_AH7LT2DSXX/SampleSheet.csv --create-fastq-for-index-reads --adapter-stringency 0 --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --output-dir /scratch/lg1883/H7LT2DSXX/bcl2fastq --fastq-compression-level 3
#     User time (seconds): 100043.79
#     System time (seconds): 984.11
#     Percent of CPU this job got: 2686%
#     Elapsed (wall clock) time (h:mm:ss or m:ss): 1:02:40
#     Average shared text size (kbytes): 0
#     Average unshared data size (kbytes): 0
#     Average stack size (kbytes): 0
#     Average total size (kbytes): 0
#     Maximum resident set size (kbytes): 22533504
#     Average resident set size (kbytes): 0
#     Major (requiring I/O) page faults: 23
#     Minor (reclaiming a frame) page faults: 152381233
#     Voluntary context switches: 821643
#     Involuntary context switches: 4676256
#     Swaps: 0
#     File system inputs: 773592064
#     File system outputs: 1043114368
#     Socket messages sent: 0
#     Socket messages received: 0
#     Signals delivered: 0
#     Page size (bytes): 4096
#     Exit status: 0

# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_deml.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_bam.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_basecall.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_bcl2fastq.sh

# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_fastq.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_cram.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_estimate.sh
