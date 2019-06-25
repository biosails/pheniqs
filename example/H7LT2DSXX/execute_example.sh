#!/usr/bin/env zsh

chmod +x H7LT2DSXX_basecall.sh
time ./H7LT2DSXX_basecall.sh

time pheniqs mux --config H7LT2DSXX_l01_estimate.json
time pheniqs mux --config H7LT2DSXX_l02_estimate.json
time pheniqs mux --config H7LT2DSXX_l03_estimate.json
time pheniqs mux --config H7LT2DSXX_l04_estimate.json
