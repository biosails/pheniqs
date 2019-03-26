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

    def update_model_prior_estimate(self):
        def update_classifier_model_prior_estimate(classifier_model, classifier_report):
            # Noise prior: {decoder low conditional confidence count} / {decoder count}
            # Barcode prior: {barcode count} + {barcode low confidence count} / {decoder count}

            classifier_estimate = {}
            classifier_model['estimate'] = classifier_estimate
            for item in [
                'count',
                'classified count',
                'classified fraction',
                'average classified confidence',
                'average classified distance',
                'pf count',
                'pf classified count',
                'pf fraction',
                'pf classified fraction',
                'classified pf fraction',
                'average pf classified confidence',
                'average pf classified distance',
                'low confidence count',
                'low conditional confidence count',
            ]:
                if item in classifier_report:
                    classifier_estimate[item] = classifier_report[item]
                else:
                    classifier_estimate[item] = 0

            # noise prior estimation
            classifier_estimate['estimated noise'] = classifier_estimate['low conditional confidence count'] / classifier_estimate['count']
            classifier_estimate['estimated noise deviation'] = 1.0 - classifier_estimate['estimated noise'] / classifier_model['simulated noise']

            barcode_report_by_hash = {}
            for barcode in classifier_report['classified']:
                hash = ''.join(barcode['barcode'])
                barcode_report_by_hash[hash] = barcode

            if 'codec' in classifier_model:
                for barcode_model in classifier_model['codec'].values():
                    hash = ''.join(barcode_model['barcode'])
                    if hash in barcode_report_by_hash:
                        barcode_report = barcode_report_by_hash[hash]
                        barcode_estimate = {}
                        barcode_model['estimate'] = barcode_estimate
                        for item in [
                            'count',
                            'average confidence',
                            'average distance',
                            'pooled classified fraction',
                            'pooled fraction',
                            'pf count',
                            'pf fraction',
                            'average pf confidence',
                            'average pf distance',
                            'pf pooled classified fraction',
                            'pf pooled fraction',
                            'low confidence count',
                        ]:
                            if item in barcode_report:
                                barcode_estimate[item] = barcode_report[item]
                            else:
                                barcode_estimate[item] = 0

                        # barcode prior estimation
                        barcode_estimate['estimated concentration'] = ((barcode_estimate['count'] + barcode_estimate['low confidence count']) / classifier_estimate['count'])
                        barcode_estimate['estimated concentration deviation'] = 1.0 - barcode_estimate['estimated concentration'] / barcode_model['simulated concentration']

            if 'unclassified' in classifier_model:
                barcode_model = classifier_model['unclassified']
                barcode_report = classifier_report['unclassified']
                barcode_estimate = {}
                for item in [
                    'count',
                    'pf count',
                    'pf fraction',
                    'pf pooled fraction',
                    'pooled classified fraction',
                    'pooled fraction',
                ]:
                    if item in barcode_report:
                        barcode_estimate[item] = barcode_report[item]
                    else:
                        barcode_estimate[item] = 0
                barcode_model['estimate'] = barcode_estimate
        def update_classifier_model_prior_estimate_new(classifier_model, classifier_report):
            classifier_estimate = {}
            classifier_model['estimate'] = classifier_estimate
            for item in [
                'count',
                'classified count',
                'classified fraction',
                'average classified confidence',
                'average classified distance',
                'pf count',
                'pf classified count',
                'pf fraction',
                'pf classified fraction',
                'classified pf fraction',
                'average pf classified confidence',
                'average pf classified distance',
                'low confidence count',
                'low conditional confidence count',
            ]:
                if item in classifier_report:
                    classifier_estimate[item] = classifier_report[item]
                else:
                    classifier_estimate[item] = 0

            # noise prior estimation
            noise_count = 0
            if 'low conditional confidence count' in classifier_estimate:
                noise_count += classifier_estimate['low conditional confidence count']

            if 'low confidence count' in classifier_estimate:
                noise_count += (classifier_estimate['low confidence count'] * (1.0 - classifier_estimate['average classified confidence']))

            classifier_estimate['estimated noise'] = noise_count / classifier_estimate['count']
            classifier_estimate['estimated noise deviation'] = 1.0 - classifier_estimate['estimated noise'] / classifier_model['simulated noise']

            self.log.info('estimating noise %s %s %s', noise_count, classifier_estimate['count'], classifier_estimate['estimated noise'])
            not_noise = 1.0 - classifier_estimate['estimated noise']
            barcode_report_by_hash = {}
            for barcode in classifier_report['classified']:
                hash = ''.join(barcode['barcode'])
                barcode_report_by_hash[hash] = barcode

            if 'codec' in classifier_model:
                for barcode_model in classifier_model['codec'].values():
                    hash = ''.join(barcode_model['barcode'])
                    if hash in barcode_report_by_hash:
                        barcode_report = barcode_report_by_hash[hash]
                        barcode_estimate = {}
                        barcode_model['estimate'] = barcode_estimate
                        for item in [
                            'count',
                            'average confidence',
                            'average distance',
                            'pooled classified fraction',
                            'pooled fraction',
                            'pf count',
                            'pf fraction',
                            'average pf confidence',
                            'average pf distance',
                            'pf pooled classified fraction',
                            'pf pooled fraction',
                            'low confidence count',
                        ]:
                            if item in barcode_report:
                                barcode_estimate[item] = barcode_report[item]
                            else:
                                barcode_estimate[item] = 0

                        # barcode prior estimation
                        barcode_estimate['estimated concentration'] = not_noise * barcode_estimate['pf pooled classified fraction']
                        barcode_estimate['estimated concentration deviation'] = 1.0 - barcode_estimate['estimated concentration'] / barcode_model['simulated concentration']

            if 'unclassified' in classifier_model:
                barcode_model = classifier_model['unclassified']
                barcode_report = classifier_report['unclassified']
                barcode_estimate = {}
                for item in [
                    'count',
                    'pf count',
                    'pf fraction',
                    'pf pooled fraction',
                    'pooled classified fraction',
                    'pooled fraction',
                ]:
                    if item in barcode_report:
                        barcode_estimate[item] = barcode_report[item]
                    else:
                        barcode_estimate[item] = 0
                barcode_model['estimate'] = barcode_estimate

        path = os.path.join(self.home, self.location['pamld prior estimate path'])
        if os.path.exists(path):
            with io.open(path, 'rb') as file:
                estimation_report = json.loads(file.read().decode('utf8'))
                for classifier_type in [ 'multiplex', 'cellular', 'molecular' ]:
                    if classifier_type in estimation_report:
                        classifier_report = estimation_report[classifier_type]
                        classifier_model = self.model[classifier_type]

                        if isinstance(classifier_model, dict):
                            update_classifier_model_prior_estimate_new(classifier_model, classifier_report)
                            self.dirty = True

                        elif isinstance(classifier_model, list):
                            model_by_index = {}
                            for item in classifier_model:
                                model_by_index[item['index']]  = item

                            for report_item in classifier_report:
                                model_item = model_by_index[report_item['index']]
                                update_classifier_model_prior_estimate_new(model_item, report_item)
                                self.dirty = True

        else: raise NoConfigurationFileError('estimation file {} not found'.format(path))

    def estimate_with_pheniqs(self):
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
            self.log.info('skipping prior estimation pheniqs execution because %s exists', self.location['pamld prior estimate path'])

    def execute(self):
        self.estimate_with_pheniqs()
        self.update_model_prior_estimate()

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
    def adjusted(self):
        return self.ontology['adjusted']

    def execute(self):
        self.instruction['output'] = os.path.join(self.home, self.location['pamld adjusted configuration path'])
        if not os.path.exists(self.instruction['output']):
            self.load_original()
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

    def adjust_prior(self):
        def adjust_decoder_prior(classifier_model, classifier_configuration):
            if 'estimate' in classifier_model:
                classifier_configuration['noise'] = classifier_model['estimate']['estimated noise']

                for barcode_key, barcode_model in classifier_model['codec'].items():
                    barcode_configuration = classifier_configuration['codec'][barcode_key]
                    if 'estimate' in barcode_model:
                        barcode_configuration['concentration'] = barcode_model['estimate']['estimated concentration']

        self.ontology['adjusted'] = deepcopy(self.original)
        for classifier_type in [ 'multiplex', 'cellular', 'molecular' ]:
            if classifier_type in self.model:
                classifier_model = self.model[classifier_type]
                if isinstance(classifier_model, dict):
                    adjust_decoder_prior(classifier_model, self.adjusted[classifier_type])

                elif isinstance(classifier_model, list):
                    model_by_index = {}
                    for item in classifier_model:
                        model_by_index[item['index']]  = item

                    for configurtion_item in self.adjusted[classifier_type]:
                        model_item = model_by_index[configurtion_item['index']]
                        update_classifier_model_prior_estimate(model_item, configurtion_item)

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
