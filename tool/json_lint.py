#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

import sys
import json

command = { 'argument': sys.argv, 'length': len(sys.argv) }

try:
    if command['length'] > 1:
        for index, argument in enumerate(command['argument']):
            if argument == '-p' or argument == '--pretty':
                print(json.dumps(json.load(sys.stdin), sort_keys=True, ensure_ascii=True, allow_nan=False,  indent=4))

            elif argument == '-c' or argument == '--compact':
                print(json.dumps(json.load(sys.stdin), sort_keys=True, ensure_ascii=True, allow_nan=False,  indent=None))
    else:
        print(json.dumps(json.load(sys.stdin), sort_keys=True, ensure_ascii=False, indent=4))

except json.decoder.JSONDecodeError as e:
    print(e)
    sys.exit(1)

sys.exit(0)
