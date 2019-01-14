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

import io
import sys
import json
from datetime import datetime

def serialize_configuration(path):
    try:
        buffer = [
            '/*',
            '    Pheniqs : PHilology ENcoder wIth Quality Statistics',
            '    Copyright (C) 2018  Lior Galanti',
            '    NYU Center for Genetics and System Biology',
            '',
            '    Author: Lior Galanti <lior.galanti@nyu.edu>',
            '',
            '    This program is free software: you can redistribute it and/or modify',
            '    it under the terms of the GNU Affero General Public License as',
            '    published by the Free Software Foundation, either version 3 of the',
            '    License, or (at your option) any later version.',
            '',
            '    This program is distributed in the hope that it will be useful,',
            '    but WITHOUT ANY WARRANTY; without even the implied warranty of',
            '    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the',
            '    GNU Affero General Public License for more details.',
            '',
            '    You should have received a copy of the GNU Affero General Public License',
            '    along with this program.  If not, see <http://www.gnu.org/licenses/>.',
            '',
            '    This file is auto generated from configuration.json',
            '',
            '    Generated on: {}'.format(datetime.now().isoformat()),
            '*/',
            '',
        ]
        buffer.append('#ifndef PHENIQS_CONFIGURATION_H')
        buffer.append('#define PHENIQS_CONFIGURATION_H')
        buffer.append('')

        with io.open(path, 'rb') as file:
            serial = json.dumps(json.loads(file.read().decode('utf8')), sort_keys=True, ensure_ascii=True, allow_nan=False,  indent=None)
            binary = [ hex(ord(c)) for c in serial ]
            length = len(binary)
            width = 0xf
            position = 0

            buffer.append('size_t configuration_json_len = {};'.format(length))
            buffer.append('')
            buffer.append('const char configuration_json[] = {')

            while position < length:
                line = '    {},'.format(', '.join(binary[ position : position + width ]))
                position += width
                if position >= length:
                    line = line[:-1]
                buffer.append(line)

            buffer.append('};')
            buffer.append('')
            buffer.append('#endif /* PHENIQS_CONFIGURATION_H */')

        print('\n'.join(buffer))

    except (OSError, json.decoder.JSONDecodeError) as e:
        print(e)
        sys.exit(1)

serialize_configuration('configuration.json')

sys.exit(0)
