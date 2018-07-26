#!/bin/bash

set -x -e

export LD_LIBRARY_PATH="${PREFIX}/lib"
make all PREFIX=${PREFIX}
make install PREFIX=${PREFIX}

#This is really a test
pheniqs demux --config test/BDGGG/BDGGG_interleave.json --validate --distance
pheniqs demux --config test/BDGGG/BDGGG_interleave.json --compile >> test/BDGGG/BDGGG.log 2>&1
pheniqs demux --config test/BDGGG/BDGGG_interleave.json >> test/BDGGG/BDGGG.log 2>&1
pheniqs demux --config test/BDGGG/BDGGG_annotated.json --validate --distance
pheniqs demux --config test/BDGGG/BDGGG_annotated.json --compile >> test/BDGGG/BDGGG.log 2>&1
pheniqs demux --config test/BDGGG/BDGGG_annotated.json >> test/BDGGG/BDGGG.log 2>&1
