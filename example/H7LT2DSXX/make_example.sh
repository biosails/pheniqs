#!/usr/bin/env zsh

pheniqs-illumina-api.py basecall illumina/181014_A00534_0024_AH7LT2DSXX
pheniqs-illumina-api.py core illumina/181014_A00534_0024_AH7LT2DSXX --no-input-npf
pheniqs-illumina-api.py multiplex illumina/181014_A00534_0024_AH7LT2DSXX
pheniqs-illumina-api.py estimate illumina/181014_A00534_0024_AH7LT2DSXX
pheniqs-illumina-api.py interleave illumina/181014_A00534_0024_AH7LT2DSXX
