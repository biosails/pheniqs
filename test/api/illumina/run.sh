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

[ -d test/api/illumina/result ] && rm -rf test/api/illumina/result;
mkdir test/api/illumina/result

illumina_basecall_test() {
  PHENIQS_TEST_NAME="illumina_basecall"
  ( cd test/api/illumina/result;
    ../../../../tool/pheniqs-illumina-api.py basecall ../181014_A00534_0024_AH7LT2DSXX;

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        PHENIQS_DIFF_FAIL="0"
        if [ "$(diff -q H7LT2DSXX_basecall_sample_sheet.csv ../valid/H7LT2DSXX_basecall_sample_sheet.csv)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected basecall sample sheet\n";
            diff H7LT2DSXX_basecall_sample_sheet.csv ../valid/H7LT2DSXX_basecall_sample_sheet.csv
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_basecall.sh ../valid/H7LT2DSXX_basecall.sh)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected basecall sample sheet\n";
            diff H7LT2DSXX_basecall.sh ../valid/H7LT2DSXX_basecall.sh
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$PHENIQS_DIFF_FAIL" != "0" ]; then
            return 1
        fi
    else
        printf "$PHENIQS_TEST_NAME returned $PHENIQS_TEST_RETURN_CODE\n";
        return $PHENIQS_TEST_RETURN_CODE
    fi
  )
};

illumina_core_test() {
  PHENIQS_TEST_NAME="illumina_core"
  ( cd test/api/illumina/result;
    ../../../../tool/pheniqs-illumina-api.py core ../181014_A00534_0024_AH7LT2DSXX;

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        PHENIQS_DIFF_FAIL="0"
        if [ "$(diff -q H7LT2DSXX_core.json ../valid/H7LT2DSXX_core.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected core configuration\n";
            diff H7LT2DSXX_core.json ../valid/H7LT2DSXX_core.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$PHENIQS_DIFF_FAIL" != "0" ]; then
            return 1
        fi
    else
        printf "$PHENIQS_TEST_NAME returned $PHENIQS_TEST_RETURN_CODE\n";
        return $PHENIQS_TEST_RETURN_CODE
    fi
  )
};

illumina_sample_test() {
  PHENIQS_TEST_NAME="illumina_sample"
  ( cd test/api/illumina/result;
    ../../../../tool/pheniqs-illumina-api.py sample ../181014_A00534_0024_AH7LT2DSXX;

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        PHENIQS_DIFF_FAIL="0"
        if [ "$(diff -q H7LT2DSXX_l01_sample.json ../valid/H7LT2DSXX_l01_sample.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected sample configuration for lane 1\n";
            diff H7LT2DSXX_l01_sample.json ../valid/H7LT2DSXX_l01_sample.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_l02_sample.json ../valid/H7LT2DSXX_l02_sample.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected sample configuration for lane 2\n";
            diff H7LT2DSXX_l02_sample.json ../valid/H7LT2DSXX_l02_sample.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_l03_sample.json ../valid/H7LT2DSXX_l03_sample.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected sample configuration for lane 3\n";
            diff H7LT2DSXX_l03_sample.json ../valid/H7LT2DSXX_l03_sample.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_l04_sample.json ../valid/H7LT2DSXX_l04_sample.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected sample configuration for lane 4\n";
            diff H7LT2DSXX_l04_sample.json ../valid/H7LT2DSXX_l04_sample.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$PHENIQS_DIFF_FAIL" != "0" ]; then
            return 1
        fi
    else
        printf "$PHENIQS_TEST_NAME returned $PHENIQS_TEST_RETURN_CODE\n";
        return $PHENIQS_TEST_RETURN_CODE
    fi
  )
};

illumina_estimate_test() {
  PHENIQS_TEST_NAME="illumina_estimate"
  ( cd test/api/illumina/result;
    ../../../../tool/pheniqs-illumina-api.py estimate ../181014_A00534_0024_AH7LT2DSXX;

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        PHENIQS_DIFF_FAIL="0"
        if [ "$(diff -q H7LT2DSXX_l01_estimate.json ../valid/H7LT2DSXX_l01_estimate.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected estimate configuration for lane 1\n";
            diff H7LT2DSXX_l01_estimate.json ../valid/H7LT2DSXX_l01_estimate.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_l02_estimate.json ../valid/H7LT2DSXX_l02_estimate.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected estimate configuration for lane 2\n";
            diff H7LT2DSXX_l02_estimate.json ../valid/H7LT2DSXX_l02_estimate.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_l03_estimate.json ../valid/H7LT2DSXX_l03_estimate.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected estimate configuration for lane 3\n";
            diff H7LT2DSXX_l03_estimate.json ../valid/H7LT2DSXX_l03_estimate.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_l04_estimate.json ../valid/H7LT2DSXX_l04_estimate.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected estimate configuration for lane 4\n";
            diff H7LT2DSXX_l04_estimate.json ../valid/H7LT2DSXX_l04_estimate.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$PHENIQS_DIFF_FAIL" != "0" ]; then
            return 1
        fi
    else
        printf "$PHENIQS_TEST_NAME returned $PHENIQS_TEST_RETURN_CODE\n";
        return $PHENIQS_TEST_RETURN_CODE
    fi
  )
};

illumina_interleave_test() {
  PHENIQS_TEST_NAME="illumina_interleave"
  ( cd test/api/illumina/result;
    ../../../../tool/pheniqs-illumina-api.py interleave ../181014_A00534_0024_AH7LT2DSXX;

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        PHENIQS_DIFF_FAIL="0"
        if [ "$(diff -q H7LT2DSXX_l01_interleave.json ../valid/H7LT2DSXX_l01_interleave.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected interleave configuration for lane 1\n";
            diff H7LT2DSXX_l01_interleave.json ../valid/H7LT2DSXX_l01_interleave.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_l02_interleave.json ../valid/H7LT2DSXX_l02_interleave.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected interleave configuration for lane 2\n";
            diff H7LT2DSXX_l02_interleave.json ../valid/H7LT2DSXX_l02_interleave.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_l03_interleave.json ../valid/H7LT2DSXX_l03_interleave.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected interleave configuration for lane 3\n";
            diff H7LT2DSXX_l03_interleave.json ../valid/H7LT2DSXX_l03_interleave.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$(diff -q H7LT2DSXX_l04_interleave.json ../valid/H7LT2DSXX_l04_interleave.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected interleave configuration for lane 4\n";
            diff H7LT2DSXX_l04_interleave.json ../valid/H7LT2DSXX_l04_interleave.json
            PHENIQS_DIFF_FAIL="1"
        fi

        if [ "$PHENIQS_DIFF_FAIL" != "0" ]; then
            return 1
        fi
    else
        printf "$PHENIQS_TEST_NAME returned $PHENIQS_TEST_RETURN_CODE\n";
        return $PHENIQS_TEST_RETURN_CODE
    fi
  )
};

illumina_basecall_test
illumina_core_test
illumina_sample_test
illumina_estimate_test
illumina_interleave_test
exit 0
