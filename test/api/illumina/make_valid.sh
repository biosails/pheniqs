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

[ -d $PHENIQS_TEST_HOME/valid ] && rm -rf $PHENIQS_TEST_HOME/valid;
mkdir $PHENIQS_TEST_HOME/valid

make_valid_illumina_basecall_test() {
    PHENIQS_TEST_NAME="illumina_basecall"
    PHENIQS_TEST_COMMAND="$PHENIQS_API_HOME/pheniqs-illumina-api.py basecall $PHENIQS_TEST_HOME/181014_A00534_0024_AH7LT2DSXX"
    PHENIQS_VALID_STDOUT="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.out"
    PHENIQS_VALID_STDERR="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.err"

    $PHENIQS_TEST_COMMAND > $PHENIQS_VALID_STDOUT 2> $PHENIQS_VALID_STDERR
    mv "$PHENIQS_TEST_HOME/181014_A00534_0024_AH7LT2DSXX/H7LT2DSXX_basecall_sample_sheet.csv" "$PHENIQS_TEST_HOME/valid/H7LT2DSXX_basecall_sample_sheet.csv"
    mv "$PHENIQS_TEST_HOME/181014_A00534_0024_AH7LT2DSXX/H7LT2DSXX_basecall.sh" "$PHENIQS_TEST_HOME/valid/H7LT2DSXX_basecall.sh"
}

make_valid_illumina_basecall_test
exit 0
