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

PHENIQS_TEST_HOME="test/api/illumina"

source "test/function.sh"

if [ ! -x "$PHENIQS_BIN" ]; then
    printf "no executable pheniqs at $PHENIQS_BIN\n";
    exit 10;
else
    $PHENIQS_BIN --version
fi

PHENIQS_TEST_basecall_sample_sheet="$PHENIQS_TEST_HOME/181014_A00534_0024_AH7LT2DSXX/H7LT2DSXX_basecall_sample_sheet.csv"
PHENIQS_TEST_basecall_script="$PHENIQS_TEST_HOME/181014_A00534_0024_AH7LT2DSXX/H7LT2DSXX_basecall.sh"

[ -f $PHENIQS_TEST_basecall_sample_sheet ] && rm -rf $PHENIQS_TEST_basecall_sample_sheet;
[ -f $PHENIQS_TEST_basecall_script ] && rm -rf $PHENIQS_TEST_basecall_script;

illumina_basecall_test() {
    PHENIQS_TEST_NAME="illumina_basecall"
    PHENIQS_TEST_COMMAND="$PHENIQS_API_HOME/pheniqs-illumina-api.py basecall $PHENIQS_TEST_HOME/181014_A00534_0024_AH7LT2DSXX"
    PHENIQS_VALID_basecall_sample_sheet="$PHENIQS_TEST_HOME/valid/H7LT2DSXX_basecall_sample_sheet.csv"
    PHENIQS_VALID_basecall_script="$PHENIQS_TEST_HOME/valid/H7LT2DSXX_basecall.sh"

    PHENIQS_TEST_STDOUT="$PHENIQS_TEST_HOME/result/$PHENIQS_TEST_NAME.out"
    PHENIQS_TEST_STDERR="$PHENIQS_TEST_HOME/result/$PHENIQS_TEST_NAME.err"
    PHENIQS_VALID_STDOUT="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.out"
    PHENIQS_VALID_STDERR="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.err"

    # execute
    $PHENIQS_TEST_COMMAND > $PHENIQS_TEST_STDOUT 2> $PHENIQS_TEST_STDERR

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        PHENIQS_DIFF_FAIL="0"
        if [ "$(diff -q $PHENIQS_TEST_basecall_sample_sheet $PHENIQS_VALID_basecall_sample_sheet)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected basecall sample sheet\n";
            diff $PHENIQS_TEST_basecall_sample_sheet $PHENIQS_VALID_basecall_sample_sheet
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q $PHENIQS_TEST_basecall_script $PHENIQS_VALID_basecall_script)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected basecall sample sheet\n";
            diff $PHENIQS_TEST_basecall_script $PHENIQS_VALID_basecall_script
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$PHENIQS_DIFF_FAIL" != "0" ]; then
            return 1
        fi
    else
        printf "$PHENIQS_TEST_NAME returned $PHENIQS_TEST_RETURN_CODE\n";
        return $PHENIQS_TEST_RETURN_CODE
    fi
}

illumina_basecall_test
exit 0
