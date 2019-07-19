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


# benchmark on prince

# Slurm Job_id=2166945 Name=bam Ended, Run time 03:14:24, COMPLETED, ExitCode 0
# real	194m15.463s
# user	3225m33.435s
# sys	1087m5.222s
#
# Slurm Job_id=2170993 Name=deml Ended, Run time 11:47:21, COMPLETED, ExitCode 0
# real	707m15.929s
# user	700m22.550s
# sys	6m21.810s
#
# Slurm Job_id=2166950 Name=basecall Ended, Run time 02:56:58, COMPLETED, ExitCode 0
#
# Slurm Job_id=2166948 Name=estimate Ended, Run time 01:25:40, COMPLETED, ExitCode 0
# real	85m36.300s
# user	808m46.322s
# sys	1388m28.835s

# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_fastq.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_cram.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_bam.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_estimate.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_deml.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_basecall.sh
# sbatch -p c01_17 --tasks-per-node=28 H7LT2DSXX_bcl2fastq.sh
