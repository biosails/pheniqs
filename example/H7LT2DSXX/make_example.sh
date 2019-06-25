#!/usr/bin/env zsh

pheniqs-illumina-api.py -v debug basecall \
--output-dir . \
--fastq-compression-level 3 \
181014_A00534_0024_AH7LT2DSXX

pheniqs-illumina-api.py -v debug core \
--base-input-url . \
--base-output-url . \
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
