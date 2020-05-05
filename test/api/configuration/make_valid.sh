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

[ -d $PHENIQS_TEST_HOME/valid ] && rm -rf $PHENIQS_TEST_HOME/valid;
mkdir $PHENIQS_TEST_HOME/valid

make_valid_configuration_header_test() {
    PHENIQS_TEST_NAME="configuration_header"
    PHENIQS_TEST_COMMAND="$PHENIQS_API_HOME/pheniqs-configuration-api.py header $PHENIQS_CONFIGURATION_PATH"
    PHENIQS_VALID_STDOUT="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.out"
    PHENIQS_VALID_STDERR="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.err"

    $PHENIQS_TEST_COMMAND > $PHENIQS_VALID_STDOUT 2> $PHENIQS_VALID_STDERR
}

make_valid_configuration_zsh_test() {
    PHENIQS_TEST_NAME="configuration_zsh"
    PHENIQS_TEST_COMMAND="$PHENIQS_API_HOME/pheniqs-configuration-api.py zsh $PHENIQS_CONFIGURATION_PATH"
    PHENIQS_VALID_STDOUT="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.out"
    PHENIQS_VALID_STDERR="$PHENIQS_TEST_HOME/valid/$PHENIQS_TEST_NAME.err"

    $PHENIQS_TEST_COMMAND > $PHENIQS_VALID_STDOUT 2> $PHENIQS_VALID_STDERR
}

make_valid_configuration_header_test
make_valid_configuration_zsh_test
exit 0
