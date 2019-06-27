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
