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
import re
import io
import json
import logging
from copy import deepcopy
from subprocess import Popen, PIPE

from core.error import *
from core import Job
from core import Shell
from core import merge
from core import to_json

class SensePrior(Shell):
    def __init__(self, ontology):
        Shell.__init__(self, ontology)
        self.log = logging.getLogger('SensePrior')

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    def execute(self):
        product = os.path.join(self.home, self.location['pamld prior estimate path'])
        if not os.path.exists(product):
            command = [ 'pheniqs', 'mux', '--sense-input', '--output', '/dev/null' ]
            command.append('--config')
            command.append(os.path.join(self.home, self.location['pamld uniform configuration path']))
            command.append('--input')
            command.append(os.path.join(self.home, self.location['simulated substitution path']))
            command.append('--report')
            command.append(os.path.join(self.home, self.location['pamld prior estimate path']))

            self.execution['command'] = ' '.join([str(i) for i in command])
            self.log.debug('executing %s', self.execution['command'])

            process = Popen(
                args=self.posix_time_command + command,
                cwd=self.current_working_directoy,
                stdout=PIPE,
                stderr=PIPE
            )
            output, error = process.communicate()
            self.execution['return code'] = process.returncode

            for line in output.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    self.execution['stdout'].append(line)

            for line in error.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    match = self.posix_time_head_ex.search(line)
                    if match:
                        for k,v in match.groupdict().items():
                            self.execution[k] = float(v)
                    else:
                        self.execution['stderr'].append(line)

            if self.execution['return code'] != 0:
                print(to_json(self.execution))
                raise CommandFailedError('pheniqs returned {} when estimating prior'.format(self.execution['return code']))
            else:
                self.dirty = True
        else:
            self.log.info('skipping prior estimation because %s exists', self.location['pamld prior estimate path'])

class AdjustPrior(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.log = logging.getLogger('AdjustPrior')

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    @property
    def original(self):
        return self.ontology['original']

    @property
    def estimate(self):
        return self.ontology['estimate']

    @property
    def adjusted(self):
        return self.ontology['adjusted']

    def execute(self):
        self.instruction['output'] = os.path.join(self.home, self.location['pamld adjusted configuration path'])
        if not os.path.exists(self.instruction['output']):
            self.load_original()
            self.load_estimate()
            self.adjust_prior()
            self.compare_to_model()
            self.save_adjusted_prior_pamld_config()
        else:
            self.log.info('skipping prior adjustment because %s exists', self.location['pamld adjusted configuration path'])

    def load_original(self):
        path = os.path.join(self.home, self.location['pamld uniform configuration path'])
        if os.path.exists(path):
            with io.open(path, 'rb') as file:
                self.ontology['original'] = json.loads(file.read().decode('utf8'))

                if 'output' in self.original:
                    del self.original['output']

                if 'feed' in self.original:
                    del self.original['feed']

                if 'multiplex' in self.original:
                    if 'output' in self.original['multiplex']:
                        del self.original['multiplex']['output']

                    if 'undetermined' in self.original['multiplex']:
                        undetermined = self.original['multiplex']['undetermined']
                        if 'feed' in undetermined:
                            del undetermined['feed']

                        if 'output' in undetermined:
                            del undetermined['output']

                    for barcode in self.original['multiplex']['codec'].values():
                        if 'feed' in barcode:
                            del barcode['feed']
                        if 'output' in barcode:
                            del barcode['output']

        else: raise NoConfigurationFileError('original configurtion file {} not found'.format(path))

    def load_estimate(self):
        path = os.path.join(self.home, self.location['pamld prior estimate path'])
        if os.path.exists(path):
            with io.open(path, 'rb') as file:
                self.ontology['estimate'] = json.loads(file.read().decode('utf8'))
        else: raise NoConfigurationFileError('estimate file {} not found'.format(path))

    def adjust_prior(self):
        def adjust_decoder_prior(estimate, decoder):
            # Noise prior: {decoder low conditional confidence count} / {decoder count}
            # Barcode prior: {barcode count} + {barcode low confidence count} / {decoder count}

            low_conditional_confidence_count = 0
            if 'low conditional confidence count' in estimate:
                low_conditional_confidence_count = estimate['low conditional confidence count']

            estimated_barcode_by_hash = {}
            decoder['noise'] = low_conditional_confidence_count / estimate['count']
            for barcode in estimate['classified']:
                hash = ''.join(barcode['barcode'])
                low_confidence_count = 0
                if 'low confidence count' in barcode:
                    low_confidence_count = barcode['low confidence count']

                barcode['prior'] = ((barcode['count'] + low_confidence_count) / estimate['count'])
                estimated_barcode_by_hash[hash] = barcode

            for barcode in decoder['codec'].values():
                hash = ''.join(barcode['barcode'])
                if hash in estimated_barcode_by_hash:
                    barcode['concentration'] = estimated_barcode_by_hash[hash]['prior']

        self.ontology['adjusted'] = deepcopy(self.original)
        for topic in [ 'multiplex', 'cellular', 'molecular' ]:
            if topic in self.estimate and topic in self.adjusted:
                if isinstance(self.adjusted[topic], dict):
                    if isinstance(self.estimate[topic], dict):
                        adjust_decoder_prior(self.estimate[topic], self.adjusted[topic])
                    else:
                        self.log.error('%s decoder structure does not match', topic)
                elif isinstance(self.adjusted[topic], list):
                    if isinstance(self.estimate[topic], list):
                        for estimate, decoder in zip(self.estimate[topic], self.adjusted[topic]):
                            adjust_decoder_prior(estimate, decoder)
                    else:
                        self.log.error('%s decoder structure does not match', topic)

    def compare_to_model(self):
        def compare_decoder_to_model(decoder, model):
            decoder['total classified prior deviation'] = 0
            decoder['simulated noise'] = model['simulated noise']
            decoder['noise estimate deviation'] = decoder['simulated noise'] - decoder['noise']

            modeled = {}
            for barcode in model['codec'].values():
                key = ''.join(barcode['barcode'])
                modeled[key] = barcode

            for barcode in decoder['codec'].values():
                key = ''.join(barcode['barcode'])
                if key in modeled:
                    barcode['simulated concentration'] = modeled[key]['simulated concentration']
                    barcode['concentration estimate deviation'] = barcode['simulated concentration'] - barcode['concentration']
                    decoder['total classified prior deviation'] += abs(barcode['concentration estimate deviation'])

            decoder['average classified prior deviation'] = decoder['total classified prior deviation'] / len(decoder['codec'])

        for topic in [ 'multiplex', 'cellular', 'molecular' ]:
            if topic in self.adjusted:
                if isinstance(self.adjusted[topic], dict):
                    compare_decoder_to_model(self.adjusted[topic], self.model[topic])

                elif isinstance(self.adjusted[topic], list):
                    for decoder, model in zip(self.adjusted[topic], self.model[topic]):
                        compare_decoder_to_model(decoder, model)

    def save_adjusted_prior_pamld_config(self):
        path = os.path.join(self.instruction['output'])
        with io.open(path, 'w') as file:
            file.write(to_json(self.adjusted))
