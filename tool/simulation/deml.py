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

import os
import io
from core.error import *
from core import merge
from core import to_json
from simulation.transcode import Transcode

class ToDeML(Transcode):
    def __init__(self, ontology):
        Transcode.__init__(self, ontology)

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    def execute(self):
        self.instruction['input'] = os.path.join(self.home, self.location['simulated substitution path'])
        self.instruction['output'] = os.path.join(self.home, self.location['deml simulated substitution path'])

        if not os.path.exists(self.instruction['output']):
            Transcode.execute(self)
        else:
            self.log.info('skipping deml substitution simulation because %s exists', self.location['deml simulated substitution path'])

    def load(self):
        Transcode.load(self)

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            XI = read['segment'][1]['fixed'][9]
            YI = read['segment'][1]['fixed'][10]
            XJ = read['segment'][2]['fixed'][9]
            YJ = read['segment'][2]['fixed'][10]
            read['segment'] = [ read['segment'][0], read['segment'][3] ]
            for segment in read['segment']:
                del segment['auxiliary']['FI']
                del segment['auxiliary']['TC']
                segment['auxiliary']['XI'] = { 'TAG': 'XI', 'TYPE': 'Z', 'VALUE': XI }
                segment['auxiliary']['YI'] = { 'TAG': 'YI', 'TYPE': 'Z', 'VALUE': YI }
                segment['auxiliary']['XJ'] = { 'TAG': 'XJ', 'TYPE': 'Z', 'VALUE': XJ }
                segment['auxiliary']['YJ'] = { 'TAG': 'YJ', 'TYPE': 'Z', 'VALUE': YJ }
            self.output_buffer.append(read)

    def summarize(self):
        pass
