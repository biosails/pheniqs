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

[ -d test/api/illumina/valid ] && rm -rf test/api/illumina/valid;
mkdir test/api/illumina/valid

make_valid_illumina_basecall_test() {
  ( cd test/api/illumina/valid;
      ../../../../tool/pheniqs-illumina-api.py basecall ../181014_A00534_0024_AH7LT2DSXX;
  )
};

make_valid_illumina_core_test() {
  ( cd test/api/illumina/valid;
      ../../../../tool/pheniqs-illumina-api.py core ../181014_A00534_0024_AH7LT2DSXX;
  )
};

make_valid_illumina_sample_test() {
  ( cd test/api/illumina/valid;
      ../../../../tool/pheniqs-illumina-api.py sample ../181014_A00534_0024_AH7LT2DSXX;
  )
};

make_valid_illumina_estimate_test() {
  ( cd test/api/illumina/valid;
      ../../../../tool/pheniqs-illumina-api.py estimate ../181014_A00534_0024_AH7LT2DSXX;
  )
};

make_valid_illumina_interleave_test() {
  ( cd test/api/illumina/valid;
      ../../../../tool/pheniqs-illumina-api.py interleave ../181014_A00534_0024_AH7LT2DSXX;
  )
};

make_valid_illumina_basecall_test
make_valid_illumina_core_test
make_valid_illumina_sample_test
make_valid_illumina_estimate_test
make_valid_illumina_interleave_test
exit 0
