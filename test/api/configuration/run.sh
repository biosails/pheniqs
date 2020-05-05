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

PHENIQS_TEST_HOME="test/api/configuration"

source "test/function.sh"

if [ ! -x "$PHENIQS_BIN" ]; then
    printf "no executable pheniqs at $PHENIQS_BIN\n";
    exit 10;
else
    $PHENIQS_BIN --version
fi

[ -d $PHENIQS_TEST_HOME/result ] && rm -rf $PHENIQS_TEST_HOME/result;
mkdir $PHENIQS_TEST_HOME/result

configuration_header_test() {
    PHENIQS_TEST_NAME="configuration_header"
    PHENIQS_TEST_COMMAND="$PHENIQS_API_HOME/pheniqs-configuration-api.py header $PHENIQS_CONFIGURATION_PATH"

    PHENIQS_TEST_STDOUT="$PHENIQS_TEST_HOME/result/$PHENIQS_TEST_NAME.out"
    PHENIQS_TEST_STDERR="$PHENIQS_TEST_HOME/result/$PHENIQS_TEST_NAME.err"
    PHENIQS_VALID_STDOUT="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.out"
    PHENIQS_VALID_STDERR="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.err"

    # execute
    $PHENIQS_TEST_COMMAND > $PHENIQS_TEST_STDOUT 2> $PHENIQS_TEST_STDERR

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        PHENIQS_DIFF_FAIL="0"
        if [ "$(diff -q -I "Generated on " $DIFF_ARGUMENTS $PHENIQS_VALID_STDOUT $PHENIQS_TEST_STDOUT)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected stdout\n";
            diff -I "Generated on " $DIFF_ARGUMENTS $PHENIQS_VALID_STDOUT $PHENIQS_TEST_STDOUT
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q $DIFF_ARGUMENTS $PHENIQS_TEST_STDERR $PHENIQS_VALID_STDERR)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected stderr\n";
            diff $DIFF_ARGUMENTS $PHENIQS_TEST_STDERR $PHENIQS_VALID_STDERR
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

configuration_zsh_test() {
    PHENIQS_TEST_NAME="configuration_zsh"
    PHENIQS_TEST_COMMAND="$PHENIQS_API_HOME/pheniqs-configuration-api.py zsh $PHENIQS_CONFIGURATION_PATH"

    PHENIQS_TEST_STDOUT="$PHENIQS_TEST_HOME/result/$PHENIQS_TEST_NAME.out"
    PHENIQS_TEST_STDERR="$PHENIQS_TEST_HOME/result/$PHENIQS_TEST_NAME.err"
    PHENIQS_VALID_STDOUT="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.out"
    PHENIQS_VALID_STDERR="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.err"

    # execute
    $PHENIQS_TEST_COMMAND > $PHENIQS_TEST_STDOUT 2> $PHENIQS_TEST_STDERR

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        PHENIQS_DIFF_FAIL="0"
        if [ "$(diff -q -I "Generated on " $DIFF_ARGUMENTS $PHENIQS_VALID_STDOUT $PHENIQS_TEST_STDOUT)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected stdout\n";
            diff -I "Generated on " $DIFF_ARGUMENTS $PHENIQS_VALID_STDOUT $PHENIQS_TEST_STDOUT
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q $DIFF_ARGUMENTS $PHENIQS_TEST_STDERR $PHENIQS_VALID_STDERR)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected stderr\n";
            diff $DIFF_ARGUMENTS $PHENIQS_TEST_STDERR $PHENIQS_VALID_STDERR
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

configuration_header_test
configuration_zsh_test
exit 0
