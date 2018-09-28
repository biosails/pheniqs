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

if [ ! -x "$PHENIQS_BIN" ]; then
    printf "no executable pheniqs at $PHENIQS_BIN\n";
    exit 10;
else
    $PHENIQS_BIN --version
fi

[ -d $PHENIQS_TEST_HOME/result ] && rm -rf $PHENIQS_TEST_HOME/result;
mkdir $PHENIQS_TEST_HOME/result

run_test \
"test/BDGGG" \
"validate_interleave" \
"demux --config test/BDGGG/BDGGG_interleave.json --precision $PHENIQS_PRECISION --validate"
PHENIQS_TEST_RETURN_CODE="$?"
if [ "$PHENIQS_TEST_RETURN_CODE" != "0" ]; then
    printf "validate_interleave failed with code $PHENIQS_TEST_RETURN_CODE\n";
    exit $PHENIQS_TEST_RETURN_CODE;
fi

run_test \
"test/BDGGG" \
"compile_interleave" \
"demux --config test/BDGGG/BDGGG_interleave.json --precision $PHENIQS_PRECISION --compile"
PHENIQS_TEST_RETURN_CODE="$?"
if [ "$PHENIQS_TEST_RETURN_CODE" != "0" ]; then
    printf "compile_interleave failed with code $PHENIQS_TEST_RETURN_CODE\n";
    exit $PHENIQS_TEST_RETURN_CODE;
fi

run_test \
"test/BDGGG" \
"validate_annotated" \
"demux --config test/BDGGG/BDGGG_annotated.json --precision $PHENIQS_PRECISION --validate --distance"
PHENIQS_TEST_RETURN_CODE="$?"
if [ "$PHENIQS_TEST_RETURN_CODE" != "0" ]; then
    printf "validate_annotated failed with code $PHENIQS_TEST_RETURN_CODE\n";
    exit $PHENIQS_TEST_RETURN_CODE;
fi

run_test \
"test/BDGGG" \
"compile_annotated" \
"demux --config test/BDGGG/BDGGG_annotated.json --precision $PHENIQS_PRECISION --compile"
PHENIQS_TEST_RETURN_CODE="$?"
if [ "$PHENIQS_TEST_RETURN_CODE" != "0" ]; then
    printf "compile_annotated failed with code $PHENIQS_TEST_RETURN_CODE\n";
    exit $PHENIQS_TEST_RETURN_CODE;
fi

run_test \
"test/BDGGG" \
"annotated" \
"demux --config test/BDGGG/BDGGG_annotated.json --precision $PHENIQS_PRECISION"
PHENIQS_TEST_RETURN_CODE="$?"
if [ "$PHENIQS_TEST_RETURN_CODE" != "0" ]; then
    printf "annotated failed with code $PHENIQS_TEST_RETURN_CODE\n";
    exit $PHENIQS_TEST_RETURN_CODE;
fi

exit 0
