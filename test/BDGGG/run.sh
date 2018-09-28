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
[ -d test/BDGGG/result ] && rm -rf test/BDGGG/result; mkdir test/BDGGG/result

run_test() {
    PHENIQS_TEST_HOME="$1"
    PHENIQS_TEST_NAME="$2"
    PHENIQS_TEST_COMMAND="$3"

    PHENIQS_TEST_STDOUT="$PHENIQS_TEST_HOME/result/$PHENIQS_TEST_NAME.out"
    PHENIQS_TEST_STDERR="$PHENIQS_TEST_HOME/result/$PHENIQS_TEST_NAME.err"
    PHENIQS_VALID_STDOUT="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.out"
    PHENIQS_VALID_STDERR="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.err"

    # execute
    $PHENIQS_BIN $PHENIQS_TEST_COMMAND > $PHENIQS_TEST_STDOUT 2> $PHENIQS_TEST_STDERR

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        if [ "$(diff -q $PHENIQS_VALID_STDOUT $PHENIQS_TEST_STDOUT)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected Pheniqs stdout\n";
            diff $PHENIQS_VALID_STDOUT $PHENIQS_TEST_STDOUT
            return 1
        fi

        if [ "$(diff -q $PHENIQS_TEST_STDERR $PHENIQS_VALID_STDERR)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected Pheniqs stderr\n";
            diff $PHENIQS_TEST_STDERR $PHENIQS_VALID_STDERR
            return 2
        fi
    else
        printf "Pheniqs returned $PHENIQS_TEST_RETURN_CODE\n";
        return $PHENIQS_TEST_RETURN_CODE
    fi
}

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
