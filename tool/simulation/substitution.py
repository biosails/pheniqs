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

# expected substitution rate : average error rate per nucleotide suggested by quality scores
# simulated substitution rate: actual simulated error rate
import os
import uuid
import math
import numpy
import logging
from copy import deepcopy

from core.error import *
from core import merge
from core import to_json
from core import decode_phred
from core import encode_phred
from core import prepare_directory
from simulation.transcode import Transcode

class SimulateSubstitution(Transcode):
    def __init__(self, ontology):
        Transcode.__init__(self, ontology)
        self.log = logging.getLogger('SimulateSubstitution')
        default = {
            'instruction': {
                'ssid': str(uuid.uuid4()),
            },
            'model': {
                'genealogy': {},
                'location': {},
            }
        }
        self.ontology = merge(default, self.ontology)
        self.genealogy['ssid'] = self.instruction['ssid']
        if 'error' in self.instruction:
            self.genealogy['error'] = self.instruction['error']

        self.location['substitution simulation home'] = os.path.join(self.location['barcode simulation home'], 'substitution', self.ssid)
        for key, relative in {
            'deml summary path' : 'deml_summary.txt',
            'deml simulated substitution path' : 'deml_simulated_substitution.bam',
            'deml demultiplex path' : 'deml_demultiplex.bam',
            'mdd demultiplex path' : 'mdd_demultiplex.bam',
            'mdd demultiplex report path' : 'mdd_demultiplex_report.json',
            'pamld adjusted configuration path' : 'pamld_adjusted_configuration.json',
            'pamld demultiplex path' : 'pamld_demultiplex.bam',
            'pamld uniform demultiplex path' : 'pamld_uniform_demultiplex.bam',
            'pamld uniform demultiplex report path' : 'pamld_uniform_demultiplex_report.json',
            'pamld demultiplex report path' : 'pamld_demultiplex_report.json',
            'pamld prior estimate path' : 'prior_estimate.json',
            'pamld accurate prior demultiplex path' : 'pamld_accurate_prior_demultiplex.bam',
            'pamld accurate prior demultiplex report path' : 'pamld_accurate_prior_demultiplex_report.json',
            'simulated substitution path' : 'simulated_substitution.bam',
        }.items():
            if key not in self.location:
                self.location[key] = os.path.join(self.substitution_home, relative)

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    @property
    def barcode_home(self):
        return self.location['barcode simulation home']

    @property
    def substitution_home(self):
        return self.location['substitution simulation home']

    def load(self):
        Transcode.load(self)

        # load the model
        for topic in [ 'multiplex', 'cellular', 'molecular' ]:
            if topic in self.model:
                if isinstance(self.model[topic], dict):
                    self.load_decoder_model(self.model[topic])
                    self.decoder_by_index.append(self.model[topic])

                elif isinstance(self.model[topic], list):
                    for decoder in self.model[topic]:
                        self.load_decoder_model(decoder)
                        self.decoder_by_index.append(decoder)

        # load accumulators
        for decoder in self.decoder_by_index:
            decoder['expected event'] = 0
            decoder['expected substitution rate compensation'] = 0
            decoder['expected substitution rate'] = 0
            decoder['unclassified']['simulated substitution'] = 0
            decoder['simulated classified substitution'] = 0
            decoder['simulated substitution'] = 0
            for barcode in decoder['codec'].values():
                barcode['simulated substitution'] = 0

        prepare_directory(os.path.join(self.home, self.location['substitution simulation home']), self.log)

    def load_decoder_model(self, decoder):
        quality_density = []
        for index, count in enumerate(decoder['quality distribution']):
            quality_density.append(count / decoder['quality distribution count'])

        decoder['barcode by index.compiled'] = [ None ] * (decoder['barcode cardinality'] + 1)
        decoder['barcode by index.compiled'][0] = decoder['unclassified']
        for barcode in decoder['codec'].values():
            decoder['barcode by index.compiled'][barcode['index']] = barcode

        if 'transform' in decoder:
            decoder['transform.compiled'] = self.parse_rule(decoder['transform'])
            if 'substitution frequency transform' in decoder:
                decoder['substitution frequency transform.compiled'] = self.parse_rule(decoder['substitution frequency transform'])

        if 'error' in self.instruction and self.instruction['error']:
            decoder['calibrated'] = True
            decoder['original expected substitution rate'] = decoder['expected substitution rate']
            decoder['requested substitution rate'] = self.instruction['error']
            phred_probability_base = pow(10.0, -0.1)

            # factor = decoder['original expected substitution rate'] / decoder['requested substitution rate']
            # offset = -10.0 * math.log10(factor)
            # for original_quality in range(64):
                # calibrated_quality = max(2, min(63, original_quality - offset)))
                # decoder['phred scale transform.compiled'].append(calibrated_quality)

            calibrated_substitution_rate = 0
            drift = 2
            normalization_factor = decoder['requested substitution rate'] / decoder['original expected substitution rate']
            while(abs(1.0 - drift) > 0.1):
                calibrated_qualit_distribution = []
                for q in range(64):
                    calibrated_qualit_distribution.append(0)

                decoder['phred scale transform.compiled'] = []
                decoder['phred scale transform.compiled'].append(0)

                self.log.info('calibrating expected error rate to %s with factor %s', decoder['requested substitution rate'], normalization_factor)

                for original_quality in range(1,64):
                    original_probability = pow(phred_probability_base, original_quality)
                    adjusted_probability = max(0, min(1, original_probability * normalization_factor))
                    calibrated_quality = max(2, min(63, round(-10.0 * math.log10(adjusted_probability))))
                    calibrated_probability = pow(phred_probability_base, calibrated_quality)

                    calibrated_qualit_distribution[calibrated_quality] += decoder['quality distribution'][original_quality]
                    calibrated_substitution_rate += decoder['quality distribution'][original_quality] * calibrated_probability

                    decoder['phred scale transform.compiled'].append(calibrated_quality)

                calibrated_substitution_rate /= decoder['simulated nucleotide']
                drift = decoder['requested substitution rate']/calibrated_substitution_rate
                normalization_factor *= drift

                self.log.info('requested substitution rate : %s', decoder['requested substitution rate'])
                self.log.info('calibrated error rate : %s', calibrated_substitution_rate)

                # for q in range(64):
                #     print(decoder['quality distribution'][q], calibrated_qualit_distribution[q])
        else:
            decoder['calibrated'] = False

        decoder['quality distribution'] = []
        for q in range(64):
            decoder['quality distribution'].append(0)

        self.load_substitution_model(decoder)

    def load_substitution_model(self, decoder):
        if 'substitution model' in decoder:
            decoder['substitution by quality.compiled'] = { 'segment': [] }
            for s in decoder['substitution model']['frequency']['segment']:
                segment = []
                for c in s['nucleotide substitution frequency']:
                    cycle = {}
                    for true, substitute in c.items():
                        nucleotide = []
                        for q in self.phred_scale:
                            correct = 1.0 - q['error']
                            quality = deepcopy(q)
                            quality['density'] = []

                            # no substitution
                            quality['density'].append({ 'nucleotide': true, 'density': correct })

                            # substitution
                            for o in substitute:
                                quality['density'].append({ 'nucleotide': o['to'], 'density': correct + o['density'] * q['error'] })

                            nucleotide.append(quality)
                        cycle[true] = nucleotide
                    segment.append(cycle)
                decoder['substitution by quality.compiled']['segment'].append(segment)
        else:
            raise BadConfigurationError('missing substitution model in decoder')

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            self.parse_qname(read)
            if len(read['segment']) != 4:
                print(to_json(read))
                exit(0)

            # convert segment SEQ and QUAL to a numpy char array
            for segment in read['segment']:
                segment['fixed'][9] = numpy.array(list(segment['fixed'][9]))
                segment['fixed'][10] = numpy.array(list(segment['fixed'][10]))

            for hint, decoder in zip(read['hint'], self.decoder_by_index):
                barcode = decoder['barcode by index.compiled'][hint]

                # simulate errors
                for token, barcode_segment_index in zip(decoder['transform.compiled']['token'], range(decoder['barcode segment cardinality'])):
                    for token_cycle, segment_cycle in zip(range(token['length']), range(token['start'], token['end'])):
                        nucleotide = read['segment'][token['segment']]['fixed'][9][segment_cycle]
                        quality = decode_phred(read['segment'][token['segment']]['fixed'][10][segment_cycle])

                        if decoder['calibrated']:
                            # calibrate the quality scores to match the expected error rate
                            quality = decoder['phred scale transform.compiled'][quality]
                            read['segment'][token['segment']]['fixed'][10][segment_cycle] = encode_phred(quality)

                        decoder['quality distribution'][quality] += 1

                        m = decoder['substitution by quality.compiled']['segment'][barcode_segment_index][token_cycle][nucleotide][quality]
                        y = m['error'] - decoder['expected substitution rate compensation']
                        t = decoder['expected substitution rate'] + y
                        compensation = (t - decoder['expected substitution rate']) - y
                        decoder['expected substitution rate'] = t

                        substitute = 'N'
                        event = numpy.random.random()
                        decoder['expected event'] += event
                        for option in m['density']:
                            if event < option['density']:
                                substitute = option['nucleotide']
                                break

                        if substitute != nucleotide:
                            barcode['simulated substitution'] += 1
                            read['segment'][token['segment']]['fixed'][9][segment_cycle] = substitute
                            # self.log.debug('event %f : substitute %s with quality %d for %s at cycle %d of segment %d', event, nucleotide, quality, substitute, token_cycle, barcode_segment_index)

            # convert segment SEQ and QUAL back to string
            for segment in read['segment']:
                segment['fixed'][9] = ''.join(segment['fixed'][9])
                segment['fixed'][10] = ''.join(segment['fixed'][10])

            self.output_buffer.append(read)

    def execute(self):
        self.instruction['input'] = os.path.join(self.home, self.location['simulated barcode path'])
        self.instruction['output'] = os.path.join(self.home, self.location['simulated substitution path'])

        if not os.path.exists(self.instruction['output']):
            Transcode.execute(self)
            self.is_model_dirty = True
        else:
            self.log.info('skipping substitution simulation because %s exists', self.location['simulated substitution path'])

    def finalize(self):
        for decoder in self.decoder_by_index:
            decoder['quality distribution count'] = 0
            decoder['quality distribution mean'] = 0
            decoder['quality distribution sd'] = 0
            for q in range(64):
                decoder['quality distribution count'] += decoder['quality distribution'][q]
                decoder['quality distribution mean'] += decoder['quality distribution'][q] * q
            decoder['quality distribution mean'] /= decoder['quality distribution count']

            for q in range(64):
                diff = q - decoder['quality distribution mean']
                decoder['quality distribution sd'] += decoder['quality distribution'][q] * diff * diff
            decoder['quality distribution sd'] /= decoder['quality distribution count']
            # del decoder['quality distribution']

            # unclassified
            unclassified = decoder['unclassified']
            unclassified['substitution rate per barcode'] = 0
            unclassified['simulated substitution rate'] = 0

            if unclassified['count'] > 0:
                unclassified['substitution rate per barcode'] = unclassified['simulated substitution'] / unclassified['count']

            if unclassified['simulated nucleotide'] > 0:
                unclassified['simulated substitution rate'] = unclassified['simulated substitution'] / unclassified['simulated nucleotide']

            # classified
            for barcode in decoder['codec'].values():
                barcode['substitution rate per barcode'] = 0
                barcode['simulated substitution rate'] = 0
                decoder['simulated classified substitution'] += barcode['simulated substitution']

                if barcode['count'] > 0:
                    barcode['substitution rate per barcode'] = barcode['simulated substitution'] / barcode['count']

                if barcode['simulated nucleotide'] > 0:
                    barcode['simulated substitution rate'] = barcode['simulated substitution'] / barcode['simulated nucleotide']

            del decoder['expected substitution rate compensation']
            decoder['simulated substitution'] = unclassified['simulated substitution'] + decoder['simulated classified substitution']

            decoder['simulated substitution rate'] = 0
            if decoder['simulated nucleotide'] > 0:
                decoder['expected event'] /= decoder['simulated nucleotide']
                decoder['expected substitution rate'] /= decoder['simulated nucleotide']
                decoder['simulated substitution rate'] = decoder['simulated substitution'] / decoder['simulated nucleotide']

            decoder['simulated classified substitution rate'] = 0
            if decoder['simulated classified nucleotide'] > 0:
                decoder['simulated classified substitution rate'] = decoder['simulated classified substitution'] / decoder['simulated classified nucleotide']

            if decoder['calibrated']:
                decoder['substitution rate deviation from expected'] = 1.0 - decoder['expected substitution rate'] / decoder['simulated substitution rate']
                decoder['substitution rate deviation from requested'] = 1.0 - decoder['requested substitution rate'] / decoder['simulated substitution rate']

                self.genealogy['multiplex calibrated'] = True
                self.genealogy['multiplex requested substitution rate'] = decoder['requested substitution rate']
                self.genealogy['multiplex expected substitution rate'] = decoder['expected substitution rate']
                self.genealogy['multiplex simulated substitution rate'] = decoder['simulated substitution rate']

                self.log.debug('substitution rate : requested %s, expected %s, simulated %s',
                    decoder['requested substitution rate'],
                    decoder['expected substitution rate'],
                    decoder['simulated substitution rate'])

                self.log.debug('substitution rate deviation : from expected %s, from requested %s',
                    decoder['substitution rate deviation from expected'],
                    decoder['substitution rate deviation from requested'])

            else:
                self.genealogy['multiplex calibrated'] = False
