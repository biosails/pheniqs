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

source "test/function.sh"

if [ ! -x "$PHENIQS_BIN" ]; then
    printf "no executable pheniqs at $PHENIQS_BIN\n";
    exit 10;
else
    $PHENIQS_BIN --version
fi

[ -d test/api/prior/result ] && rm -rf test/api/prior/result;
mkdir test/api/prior/result

prior_test() {
  PHENIQS_TEST_NAME="prior"
  ( cd test/api/prior/result;
    ../../../../tool/pheniqs-prior-api.py --configuration ../BDGGG_annotated.json --report ../BDGGG_annotated_report.json > BDGGG_annotated_estimated.json;

    PHENIQS_TEST_RETURN_CODE=$?
    if [ "$PHENIQS_TEST_RETURN_CODE" == "0" ]; then
        PHENIQS_DIFF_FAIL="0"
        if [ "$(diff -q BDGGG_annotated_estimated.json ../valid/BDGGG_annotated_estimated.json)" ]; then
            printf "$PHENIQS_TEST_NAME : Unexpected prior configuration\n";
            diff BDGGG_annotated_estimated.json ../valid/BDGGG_annotated_estimated.json
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

prior_test
exit 0
