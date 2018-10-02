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

source "test/BDGGG/function.sh"

[ -d $PHENIQS_TEST_HOME/valid ] && rm -rf $PHENIQS_TEST_HOME/valid;
mkdir $PHENIQS_TEST_HOME/valid

make_valid_test $PHENIQS_TEST_HOME "validate_interleave" \
"demux --config $PHENIQS_TEST_HOME/BDGGG_interleave.json --precision $PHENIQS_PRECISION --validate"

make_valid_test $PHENIQS_TEST_HOME "compile_interleave" \
"demux --config $PHENIQS_TEST_HOME/BDGGG_interleave.json --precision $PHENIQS_PRECISION --compile"

make_valid_test $PHENIQS_TEST_HOME "validate_annotated" \
"demux --config $PHENIQS_TEST_HOME/BDGGG_annotated.json --precision $PHENIQS_PRECISION --validate --distance"

make_valid_test $PHENIQS_TEST_HOME "compile_annotated" \
"demux --config $PHENIQS_TEST_HOME/BDGGG_annotated.json --precision $PHENIQS_PRECISION --compile"

make_valid_test $PHENIQS_TEST_HOME "annotated" \
"demux --config $PHENIQS_TEST_HOME/BDGGG_annotated.json --precision $PHENIQS_PRECISION"
