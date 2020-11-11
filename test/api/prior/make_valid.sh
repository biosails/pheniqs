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

[ -d test/api/prior/valid ] && rm -rf test/api/prior/valid;
mkdir test/api/prior/valid

make_valid_prior_test() {
  ( cd test/api/prior/valid;
    ../../../../tool/pheniqs-prior-api.py --configuration ../BDGGG_annotated.json --report ../BDGGG_annotated_report.json > BDGGG_annotated_estimated.json
  )
};

make_valid_prior_test
exit 0
