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

from core import *

class Prior(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)

    @property
    def original(self):
        return self.ontology['original']

    @property
    def report(self):
        return self.ontology['report']

    @property
    def adjusted(self):
        return self.ontology['adjusted']

    @property
    def model(self):
        return self.ontology['model']

    def load_model(self):
        if 'model' in self.instruction and self.instruction['model']:
            if os.path.exists(self.instruction['model']):
                self.log.debug('loading model %s', self.instruction['model'])
                with io.open(self.instruction['model'], 'rb') as file:
                    self.ontology['model'] = json.loads(file.read().decode('utf8'))
                    return True
        return False

    def load_original(self):
        if os.path.exists(self.instruction['configuration']):
            self.log.debug('loading original configuration %s', self.instruction['configuration'])
            with io.open(self.instruction['configuration'], 'rb') as file:
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

        else: raise NoConfigurationFileError('original configurtion file {} not found'.format(self.instruction['configuration']))

    def load_report(self):
        if os.path.exists(self.instruction['report']):
            self.log.debug('loading report %s', self.instruction['report'])
            with io.open(self.instruction['report'], 'rb') as file:
                self.ontology['report'] = json.loads(file.read().decode('utf8'))
        else: raise NoConfigurationFileError('report file {} not found'.format(self.instruction['report']))

    def save_adjusted(self):
        self.log.debug('saving adjusted configurtion %s', self.instruction['adjusted'])
        with io.open(self.instruction['adjusted'], 'w') as file:
            file.write(to_json(self.adjusted))

    def adjust_decoder_prior(self, report, decoder):
        reported_barcode_by_hash = {}
        decoder['noise'] = report['low conditional confidence count'] / report['count']
        for barcode in report['classified']:
            hash = ''.join(barcode['barcode'])
            low_confidence_count = 0
            if 'low confidence count' in barcode:
                low_confidence_count = barcode['low confidence count']
            barcode['prior'] = ((barcode['count'] + low_confidence_count) / report['count'])
            reported_barcode_by_hash[hash] = barcode

        for barcode in decoder['codec'].values():
            hash = ''.join(barcode['barcode'])
            if hash in reported_barcode_by_hash:
                barcode['concentration'] = reported_barcode_by_hash[hash]['prior']

    def adjust_prior(self):
        self.ontology['adjusted'] = deepcopy(self.original)
        for topic in [ 'multiplex', 'cellular', 'molecular' ]:
            if topic in self.report and topic in self.adjusted:
                if isinstance(self.adjusted[topic], dict):
                    if isinstance(self.report[topic], dict):
                        self.adjust_decoder_prior(self.report[topic], self.adjusted[topic])
                    else:
                        self.log.error('%s decoder structure does not match', topic)
                elif isinstance(self.adjusted[topic], list):
                    if isinstance(self.report[topic], list):
                        for report, decoder in zip(self.report[topic], self.adjusted[topic]):
                            self.adjust_decoder_prior(report, decoder)
                    else:
                        self.log.error('%s decoder structure does not match', topic)

    def compare_decoder_to_model(self, decoder, model):
        noise_prior_deviation = 0
        total_classified_prior_deviation = 0
        decoder['noise deviation'] = model['simulated noise'] - decoder['noise']
        noise_prior_deviation = decoder['noise deviation']

        modeled = {}
        for barcode in model['codec'].values():
            key = ''.join(barcode['barcode'])
            modeled[key] = barcode

        for barcode in decoder['codec'].values():
            key = ''.join(barcode['barcode'])
            if key in modeled:
                barcode['concentration estimate deviation'] = modeled[key]['simulated concentration'] - barcode['concentration']
                total_classified_prior_deviation += barcode['concentration estimate deviation']

        print('noise prior deviation : {}'.format(noise_prior_deviation))
        print('total classified prior deviation : {}'.format(total_classified_prior_deviation))
        print('average classified prior deviation : {}'.format(total_classified_prior_deviation / len(decoder['codec'])))

    def compare_to_model(self, decoder_key):
        if self.load_model():
            for topic in [ 'multiplex', 'cellular', 'molecular' ]:
                if topic in self.adjusted:
                    if isinstance(self.adjusted[topic], dict):
                        self.compare_decoder_to_model(self.adjusted[topic], self.model[topic])

                    elif isinstance(self.adjusted[topic], list):
                        for decoder, model in zip(self.adjusted[topic], self.model[topic]):
                            self.compare_decoder_to_model(decoder, model)

    def execute(self):
        self.load_original()
        self.load_report()
        self.adjust_prior()
        self.compare_to_model('A5KVK')
        self.save_adjusted()
