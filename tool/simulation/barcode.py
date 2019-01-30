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
import uuid
import json
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

class SimulateBarcode(Transcode):
    def __init__(self, ontology):
        Transcode.__init__(self, ontology)
        self.log = logging.getLogger('SimulateBarcode')
        default = {
            'instruction': {
                'bsid': str(uuid.uuid4()),
            },
            'vocabulary': {
                'decoder': {},
                'preset': {},
                'substitution model': {}
            },
            'model': {
                'genealogy': {},
                'location': {},
                'run count': 0,
                'noise model': {},
            },
            'execution': {},
        }
        self.ontology = merge(default, self.ontology)
        self.genealogy['bsid'] = self.instruction['bsid']
        self.location['barcode simulation home'] = os.path.join('simulation', self.bsid)
        self.location['simulated barcode path'] = os.path.join(self.location['barcode simulation home'], 'simulated_barcode.bam')
        self.location['pamld uniform configuration path'] = os.path.join(self.location['barcode simulation home'], 'pamld_uniform_configuration.json')
        self.location['mdd configuration path'] = os.path.join(self.location['barcode simulation home'], 'mdd_configuration.json')
        self.location['pamld accurate prior configuration path'] = os.path.join(self.location['barcode simulation home'], 'pamld_accurate_prior_configuration.json')
        self.location['deml index path'] = os.path.join(self.location['barcode simulation home'], 'deml_index.txt')

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def vocabulary(self):
        return self.ontology['vocabulary']

    def load(self):
        Transcode.load(self)
        self.load_vocabulary()
        self.load_preset()
        self.load_noise_model()
        prepare_directory(os.path.join(self.home, self.location['barcode simulation home']), self.log)

    def load_vocabulary(self):
        path = os.path.realpath(os.path.join(os.path.dirname(__file__), '../configuration/vocabulary.json'))
        if os.path.exists(path):
            self.log.debug('loading vocabulary %s', path)
            with io.open(path, 'rb') as file:
                self.ontology['vocabulary'] = json.loads(file.read().decode('utf8'))
        else: raise NoConfigurationFileError('vocabulary file {} not found'.format(path))

    def load_preset(self):
        def load_decoder_preset(decoder):
            compiled = decoder
            if 'base' in decoder:
                if 'decoder' in self.vocabulary and decoder['base'] in self.vocabulary['decoder']:
                    compiled = merge(self.vocabulary['decoder'][decoder['base']], decoder)
                else:
                    raise BadConfigurationError('reference to unknown base decoder {}'.format(decoder['base']))

            if 'unclassified' not in compiled: compiled['unclassified'] = {}
            compiled['unclassified']['index'] = 0
            compiled['unclassified']['description'] = 'unclassified'
            compiled['unclassified']['configured concentration'] = compiled['configured noise']
            compiled['barcode length'] = []
            compiled['nucleotide cardinality'] = 0
            compiled['barcode cardinality'] = len(compiled['codec'])

            if 'transform' in compiled:
                compiled['transform.compiled'] = self.parse_rule(compiled['transform'])
                for token in compiled['transform.compiled']['token']:
                    compiled['barcode length'].append(token['length'])
                    compiled['nucleotide cardinality'] += token['length']
                compiled['barcode segment cardinality'] = len(compiled['barcode length'])

                if 'substitution frequency transform' in compiled:
                    compiled['substitution frequency transform.compiled'] = self.parse_rule(compiled['substitution frequency transform'])
                    if compiled['barcode length'] != [ t['length'] for t in compiled['substitution frequency transform.compiled']['token'] ]:
                        raise BadConfigurationError('substitution frequency transform token array size does not match decoder {} token array'.format(str(compiled['index'])))
            else:
                raise BadConfigurationError('missing transform for decoder {}'.format(str(compiled['index'])))

            barcode_index = 1
            total_concentration = 0
            not_noise = 1.0 - compiled['configured noise']
            for barcode_key in sorted(compiled['codec'].keys()):
                barcode = compiled['codec'][barcode_key]
                barcode['key'] = barcode_key
                barcode['index'] = barcode_index
                barcode['description'] = '{}:{}'.format(barcode_index, '-'.join(barcode['barcode']))
                if 'configured concentration' not in barcode: barcode['configured concentration'] = 1.0
                total_concentration += barcode['configured concentration']

                barcode_length = [ len(i) for i in barcode['barcode'] ]
                if barcode_length != compiled['barcode length']:
                    raise BadConfigurationError('found barcode length {} for barcode {} in decoder {} and expected {}'.format(
                        ' '.join([str(i) for i in barcode_length]),
                        barcode_key,
                        str(compiled['index']),
                        ' '.join([str(i) for i in compiled['barcode length']])
                    ))

                # convert barcode sequence to a numpy char array
                barcode['barcode.compiled'] = []
                for segment in barcode['barcode']:
                    barcode['barcode.compiled'].append(numpy.array(list(segment)))
                barcode_index += 1

            density = 0
            compiled['barcode distribution.compiled'] = []
            for barcode_key, barcode_value in compiled['codec'].items():
                barcode_value['configured concentration'] = not_noise * (barcode_value['configured concentration'] / total_concentration)
                density += barcode_value['configured concentration']
                compiled['barcode distribution.compiled'].append({'density': density, 'item': barcode_value})

            if compiled['substitution model reference'] in self.vocabulary['substitution model']:
                compiled['substitution model'] = { 'frequency': { 'segment': [] } }
                for token_index, token in enumerate(compiled['substitution frequency transform.compiled']['token']):
                    compiled['substitution model']['frequency']['segment'].append({'nucleotide substitution frequency': []})

                reference = self.vocabulary['substitution model'][compiled['substitution model reference']]
                for token_index, token in enumerate(compiled['substitution frequency transform.compiled']['token']):
                    for cycle in range(token['start'], token['end']):
                        compiled['substitution model']['frequency']['segment'][token_index]['nucleotide substitution frequency'].append(deepcopy(reference['segment'][token['segment']]['nucleotide substitution frequency'][cycle]))
            else:
                raise BadConfigurationError('unknown substitution model {}'.format(decoder['substitution model reference']))

            compiled['quality distribution'] = []
            for q in range(64):
                compiled['quality distribution'].append(0)

            compiled['quality distribution by base'] = {}
            for base in ['A', 'T', 'C', 'G', 'N']:
                compiled['quality distribution by base'][base] = []
                for q in range(64):
                    compiled['quality distribution by base'][base].append(0)

            return compiled

        if self.instruction['preset'] in self.vocabulary['preset']:
            self.ontology['model'] = merge(self.ontology['model'], self.vocabulary['preset'][self.instruction['preset']])

            decoder_index = 0
            for topic in [ 'multiplex', 'cellular', 'molecular' ]:
                if topic in self.model:
                    if isinstance(self.model[topic], dict):
                        self.model[topic]['index'] = decoder_index
                        self.model[topic] = load_decoder_preset(self.model[topic])
                        self.decoder_by_index.append(self.model[topic])
                        decoder_index +=1

                    elif isinstance(self.model[topic], list):
                        buffer = self.model[topic]
                        self.model[topic] = []
                        for decoder in buffer:
                            decoder['index'] = decoder_index
                            compiled = load_decoder_preset(decoder)
                            self.model[topic].append(compiled)
                            self.decoder_by_index.append(compiled)
                            decoder_index +=1

            # load accumulators
            self.model['count'] = 0
            for decoder in self.decoder_by_index:
                decoder['classified count'] = 0
                decoder['unclassified']['count'] = 0
                decoder['simulated nucleotide'] = 0
                decoder['expected substitution rate compensation'] = 0
                decoder['expected substitution rate'] = 0
                for barcode in decoder['codec'].values():
                    barcode['count'] = 0

        else: raise BadConfigurationError('unknown preset %s', self.vocabulary['preset'])

    def load_noise_model(self):
        if 'noise model' in self.model:
            if 'reference fasta path' in self.model['noise model']:
                path = os.path.realpath(os.path.join(os.path.dirname(__file__), '../configuration', self.model['noise model']['reference fasta path']))
                if os.path.exists(path):
                    with io.open(path, 'rb') as file:
                        content = file.read().decode('ascii').split('\n')
                        if content[0][0] == '>':
                            title = content[0][1:].split(' ')
                            self.model['noise model']['genome id'] = title[0]
                            self.model['noise model']['genome comment'] = ' '.join(title[1:])
                            self.model['noise model']['genome.compiled'] = numpy.array(list(''.join([line.strip() for line in content[1:]])))
                            self.model['noise model']['nucleotide cardinality'] = len(self.model['noise model']['genome.compiled'])
                else: raise NoConfigurationFileError('noise model reference fasta {} not found'.format(self.model['noise model']['reference fasta path']))
            else: raise BadConfigurationError('missing reference fasta path for noise model')

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            read['index'] = self.increment_read_counter()
            read['hint'] = []

            # convert segment SEQ to a numpy char array
            for segment in read['segment']:
                segment['fixed'][9] = numpy.array(list(segment['fixed'][9]))

            for decoder in self.decoder_by_index:
                barcode = self.pick_barcode(decoder)
                barcode['count'] += 1
                read['hint'].append(barcode['index'])

                if barcode['index'] > 0:
                    # barcode: copy the selected barcode segment into the read
                    for token, barcode_segment in zip(decoder['transform.compiled']['token'], barcode['barcode.compiled']):
                        read['segment'][token['segment']]['fixed'][9][token['start']:token['end']] = barcode_segment

                else:
                    # noise: insert a random piece of PhiX
                    for token, segment_length in zip(decoder['transform.compiled']['token'], decoder['barcode length']):
                        offset = self.sample_random_noise(segment_length)
                        read['segment'][token['segment']]['fixed'][9][token['start']:token['end']] = self.model['noise model']['genome.compiled'][offset:offset + segment_length]

                for token, barcode_segment_index in zip(decoder['transform.compiled']['token'], range(decoder['barcode segment cardinality'])):
                    for token_cycle, segment_cycle in zip(range(token['length']), range(token['start'], token['end'])):
                        nucleotide = read['segment'][token['segment']]['fixed'][9][segment_cycle]
                        quality = decode_phred(read['segment'][token['segment']]['fixed'][10][segment_cycle])
                        decoder['quality distribution'][quality] += 1
                        decoder['quality distribution by base'][nucleotide][quality] += 1

                        y = self.phred_scale[quality]['error'] - decoder['expected substitution rate compensation']
                        t = decoder['expected substitution rate'] + y
                        compensation = (t - decoder['expected substitution rate']) - y
                        decoder['expected substitution rate'] = t

            self.simulate_qname(read)

            # convert segment SEQ back to string
            for segment in read['segment']:
                segment['fixed'][9] = ''.join(segment['fixed'][9])

            self.output_buffer.append(read)

    def increment_read_counter(self):
        self.model['count'] += 1
        return self.model['count']

    def pick_barcode(self, decoder):
        barcode = None
        event = numpy.random.random()
        for option in decoder['barcode distribution.compiled']:
            if event < option['density']:
                barcode = option['item']
                # self.log.debug('barcode %s picked by %.15f', barcode['description'], event)
                break
        if barcode is None:
            barcode = decoder['unclassified']
            # self.log.debug('unclassified picked %.15f', event)
        return barcode

    def sample_random_noise(self, size):
        scale = self.model['noise model']['nucleotide cardinality'] - size
        event = numpy.random.random()
        return int(event * scale)

    def simulate_qname(self, read):
        # @PHENIQS:1:HABDFADXX:0112241932:23:31:4
        #<instrument id>:<run count>:<flowcell id>:<read index>[:<decoder hint>]*
        qname = [ self.instrument_id, str(self.run_count), self.flowcell_id ]
        qname.append(str(read['index']).rjust(10,'0'))
        for hint in read['hint']:
            qname.append(str(hint))

        simulated_qname = ':'.join(qname)
        for segment in read['segment']:
            segment['fixed'][0] = simulated_qname

    def execute(self):
        self.instruction['input'] = os.path.realpath(os.path.expanduser(os.path.expandvars(self.instruction['input'])))
        self.instruction['output'] = os.path.join(self.home, self.location['simulated barcode path'])

        if not os.path.exists(self.instruction['output']):
            Transcode.execute(self)
            self.ontology['persistence']['dirty'] = True
        else:
            self.log.info('skipping barcode simulation because %s exists', self.location['simulated barcode path'])

        self.save_uniform_pamld_config()
        self.save_mdd_config()
        self.save_accurate_prior_pamld_config()
        self.save_deml_index()

    def finalize(self):
        for decoder in self.decoder_by_index:
            del decoder['expected substitution rate compensation']

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

            # unclassified
            unclassified = decoder['unclassified']
            unclassified['simulated concentration'] = 0
            unclassified['simulated concentration deviation'] = 0
            unclassified['simulated nucleotide'] = unclassified['count'] * decoder['nucleotide cardinality']
            if self.model['count'] > 0:
                unclassified['simulated concentration'] = unclassified['count'] / self.model['count']
            if unclassified['configured concentration'] > 0:
                unclassified['simulated concentration deviation'] = 1.0 - unclassified['simulated concentration'] / unclassified['configured concentration']

            # classified
            for barcode in decoder['codec'].values():
                barcode['simulated concentration'] = 0
                barcode['simulated concentration deviation'] = 0

                decoder['classified count'] += barcode['count']
                barcode['simulated nucleotide'] = barcode['count'] * decoder['nucleotide cardinality']
                if self.model['count'] > 0:
                    barcode['simulated concentration'] = barcode['count'] / self.model['count']
                if barcode['configured concentration'] > 0:
                    barcode['simulated concentration deviation'] = 1.0 - barcode['simulated concentration'] / barcode['configured concentration']

            decoder['simulated noise'] = 0
            decoder['simulated noise deviation'] = 0
            decoder['simulated classified nucleotide'] = decoder['classified count'] * decoder['nucleotide cardinality']
            decoder['simulated nucleotide'] = unclassified['simulated nucleotide'] + decoder['simulated classified nucleotide']
            if decoder['simulated nucleotide'] > 0:
                decoder['expected substitution rate'] /= decoder['simulated nucleotide']

            if self.model['count'] > 0:
                decoder['simulated noise'] = unclassified['count'] / self.model['count']
            if decoder['configured noise'] > 0:
                decoder['simulated noise deviation'] = 1.0 - decoder['simulated noise'] / decoder['configured noise']

    def save_deml_index(self):
        path = os.path.join(self.home, self.location['deml index path'])
        if not os.path.exists(path):
            if 'multiplex' in self.model:
                content = []
                content.append('#{}'.format('\t'.join(['Index1', 'Index2', 'Name'])))
                for k,v in self.model['multiplex']['codec'].items():
                    line = [ str(s) for s in v['barcode'] ]
                    name = ':'.join([ self.model['flowcell id'], ''.join(line) ])
                    line.append(name)
                    content.append('\t'.join(line))

                self.log.info('saving deml index %s', self.location['deml index path'])
                with io.open(path, 'w') as file:
                    file.write('\n'.join(content))
        else:
            self.log.info('skipping deml index because %s exists', self.location['deml index path'])

    def save_accurate_prior_pamld_config(self):
        def extract_codec(decoder):
            compiled = {
                'noise': decoder['simulated noise'],
                'transform': decoder['transform'],
                'algorithm': 'pamld',
                'codec': {},
            }
            for k in [
                'confidence threshold'
            ]:
                if k in decoder:
                    compiled[k] = decoder[k]

            for k,v in decoder['codec'].items():
                compiled['codec'][k] = {
                    'barcode': v['barcode'],
                    'concentration': v['simulated concentration'],
                }
            return compiled

        path = os.path.join(self.home, self.location['pamld accurate prior configuration path'])
        if not os.path.exists(path):
            ontology = {
                'PL': self.model['PL'],
                'flowcell id': self.model['flowcell id'],
                'transform': self.model['transform'],
                'include filtered': True,
            }
            for topic in [ 'multiplex', 'cellular', 'molecular' ]:
                if topic in self.model:
                    if isinstance(self.model[topic], dict):
                        ontology[topic] = extract_codec(self.model[topic])

                    elif isinstance(self.model[topic], list):
                        ontology[topic] = []
                        for decoder in self.model[topic]:
                            ontology[topic].append(extract_codec(decoder))

            self.log.info('saving accurate prior pamld configuration %s', self.location['pamld accurate prior configuration path'])
            with io.open(path, 'w') as file:
                file.write(to_json(ontology))
        else:
            self.log.info('skipping accurate prior configuration because %s exists', self.location['pamld accurate prior configuration path'])

    def save_uniform_pamld_config(self):
        def extract_codec(decoder):
            compiled = {
                'transform': decoder['transform'],
                'algorithm': 'pamld',
                'codec': {}
            }
            for k in [
                'confidence threshold'
            ]:
                if k in decoder:
                    compiled[k] = decoder[k]

            for k,v in decoder['codec'].items():
                compiled['codec'][k] = { 'barcode': v['barcode'] }
            return compiled

        path = os.path.join(self.home, self.location['pamld uniform configuration path'])
        if not os.path.exists(path):
            ontology = {
                'PL': self.model['PL'],
                'flowcell id': self.model['flowcell id'],
                'transform': self.model['transform'],
                'include filtered': True,
            }
            for topic in [ 'multiplex', 'cellular', 'molecular' ]:
                if topic in self.model:
                    if isinstance(self.model[topic], dict):
                        ontology[topic] = extract_codec(self.model[topic])

                    elif isinstance(self.model[topic], list):
                        ontology[topic] = []
                        for decoder in self.model[topic]:
                            ontology[topic].append(extract_codec(decoder))

            self.log.info('saving uniform pamld configuration %s', self.location['pamld uniform configuration path'])
            with io.open(path, 'w') as file:
                file.write(to_json(ontology))
        else:
            self.log.info('skipping uniform configuration because %s exists', self.location['pamld uniform configuration path'])

    def save_mdd_config(self):
        def extract_codec(decoder):
            compiled = {
                'transform': decoder['transform'],
                'algorithm': 'mdd',
                'codec': {}
            }
            for k,v in decoder['codec'].items():
                compiled['codec'][k] = { 'barcode': v['barcode'] }
            return compiled

        path = os.path.join(self.home, self.location['mdd configuration path'])
        if not os.path.exists(path):
            ontology = {
                'PL': self.model['PL'],
                'flowcell id': self.model['flowcell id'],
                'transform': self.model['transform'],
                'include filtered': True,
            }
            for topic in [ 'multiplex', 'cellular', 'molecular' ]:
                if topic in self.model:
                    if isinstance(self.model[topic], dict):
                        ontology[topic] = extract_codec(self.model[topic])

                    elif isinstance(self.model[topic], list):
                        ontology[topic] = []
                        for decoder in self.model[topic]:
                            ontology[topic].append(extract_codec(decoder))

            self.log.info('saving mdd configuration %s', self.location['mdd configuration path'])
            with io.open(path, 'w') as file:
                file.write(to_json(ontology))
        else:
            self.log.info('skipping mdd configuration because %s exists', self.location['mdd configuration path'])
