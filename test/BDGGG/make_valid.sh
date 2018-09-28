#!/usr/bin/env bash

# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2018  Lior Galanti
# NYU Center for Genetics and System Biology

# Author: Lior Galanti <lior.galanti@nyu.edu>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

PHENIQS_PRECISION="10"
PHENIQS_BIN="./pheniqs"
PHENIQS_TEST_HOME="test/BDGGG"

[ -d $PHENIQS_TEST_HOME/valid ] && rm -rf $PHENIQS_TEST_HOME/valid;
mkdir $PHENIQS_TEST_HOME/valid

$PHENIQS_BIN demux --config $PHENIQS_TEST_HOME/BDGGG_interleave.json --precision $PHENIQS_PRECISION --validate             > $PHENIQS_TEST_HOME/valid/validate_interleave.out 2> $PHENIQS_TEST_HOME/valid/validate_interleave.err
$PHENIQS_BIN demux --config $PHENIQS_TEST_HOME/BDGGG_interleave.json --precision $PHENIQS_PRECISION --compile              > $PHENIQS_TEST_HOME/valid/compile_interleave.out  2> $PHENIQS_TEST_HOME/valid/compile_interleave.err
$PHENIQS_BIN demux --config $PHENIQS_TEST_HOME/BDGGG_annotated.json  --precision $PHENIQS_PRECISION --validate --distance  > $PHENIQS_TEST_HOME/valid/validate_annotated.out  2> $PHENIQS_TEST_HOME/valid/validate_annotated.err
$PHENIQS_BIN demux --config $PHENIQS_TEST_HOME/BDGGG_annotated.json  --precision $PHENIQS_PRECISION --compile              > $PHENIQS_TEST_HOME/valid/compile_annotated.out   2> $PHENIQS_TEST_HOME/valid/compile_annotated.err
$PHENIQS_BIN demux --config $PHENIQS_TEST_HOME/BDGGG_annotated.json  --precision $PHENIQS_PRECISION                        > $PHENIQS_TEST_HOME/valid/annotated.out           2> $PHENIQS_TEST_HOME/valid/annotated.err
