#!/usr/bin/env zsh

chmod +x H7LT2DSXX_basecall.sh
time ./H7LT2DSXX_basecall.sh

time pheniqs mux --config H7LT2DSXX_l01_estimate.json
time pheniqs mux --config H7LT2DSXX_l02_estimate.json
time pheniqs mux --config H7LT2DSXX_l03_estimate.json
time pheniqs mux --config H7LT2DSXX_l04_estimate.json

# pheniqs mux --config H7LT2DSXX_l01_estimate.json  74530.93s user 25033.40s system 1556% cpu 1:46:37.06 total
# pheniqs mux --config H7LT2DSXX_l02_estimate.json  73217.38s user 25545.95s system 1598% cpu 1:43:00.30 total
# pheniqs mux --config H7LT2DSXX_l03_estimate.json  75480.70s user 25363.05s system 1570% cpu 1:46:59.32 total
