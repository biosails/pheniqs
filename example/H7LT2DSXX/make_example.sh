#!/usr/bin/env zsh

pheniqs-illumina-api.py -v debug basecall \
--output-dir . \
--fastq-compression-level 3 \
181014_A00534_0024_AH7LT2DSXX

pheniqs-illumina-api.py -v debug core \
--base-input . \
--base-output . \
--no-input-npf \
181014_A00534_0024_AH7LT2DSXX

pheniqs-illumina-api.py -v debug multiplex \
--confidence 0.95 \
--noise 0.05 \
181014_A00534_0024_AH7LT2DSXX

pheniqs-illumina-api.py -v debug estimate \
--confidence 0.95 \
--noise 0.05 \
181014_A00534_0024_AH7LT2DSXX

pheniqs-illumina-api.py -v debug interleave \
181014_A00534_0024_AH7LT2DSXX

pheniqs-prior-api.py \
--report H7LT2DSXX_l01_estimate_report.json \
--configuration H7LT2DSXX_l01_sample.json \
> H7LT2DSXX_l01_adjusted.json

pheniqs-prior-api.py \
--report H7LT2DSXX_l02_estimate_report.json \
--configuration H7LT2DSXX_l02_sample.json \
> H7LT2DSXX_l02_adjusted.json

pheniqs-prior-api.py \
--report H7LT2DSXX_l03_estimate_report.json \
--configuration H7LT2DSXX_l03_sample.json \
> H7LT2DSXX_l03_adjusted.json

pheniqs-prior-api.py \
--report H7LT2DSXX_l04_estimate_report.json \
--configuration H7LT2DSXX_l04_sample.json \
> H7LT2DSXX_l04_adjusted.json
