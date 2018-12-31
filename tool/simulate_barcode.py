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

from transcode import *

class SimulateBarcode(TranscodeSAM):
    def __init__(self, ontology):
        TranscodeSAM.__init__(self, ontology)

    @property
    def configuration(self):
        return self.ontology['configuration']

    @property
    def noise_model(self):
        return self.model['noise model']

    @property
    def noise_genome(self):
        return self.model['noise model']['genome.compiled']

    def report(self):
        if 'report' in self.instruction and self.instruction['report']:
            report = deepcopy(self.model)
            remove_compiled(report)
            report['genealogy']['document type'] = 'barcode simulation report';
            with io.open(self.instruction['report'], 'w') as file:
                file.write(to_json(report))

    def load(self):
        self.load_configuration()
        self.load_preset()
        self.load_noise_model()
        self.load_report()

    def load_configuration(self):
        if os.path.exists(self.instruction['configuration']):
            self.instruction['working directory'] = os.path.dirname(self.instruction['configuration'])
            self.log.debug('loading configuration %s', self.instruction['configuration'])
            with io.open(self.instruction['configuration'], 'rb') as file:
                self.ontology['configuration'] = json.loads(file.read().decode('utf8'))
        else: raise NoConfigurationFileError('configure {} not found'.format(self.instruction['configuration']))

    def load_preset(self):
        if self.instruction['preset'] in self.configuration['preset']:
            self.ontology['model'] = {
                'run count': 0,
                'genealogy': {
                    'document type': 'barcode simulation model',
                    'barcode simulation id': str(uuid.uuid4()),
                }
            }
            self.ontology['model'] = merge(self.ontology['model'], self.configuration['preset'][self.instruction['preset']])

            decoder_index = 0
            for topic in [ 'multiplex', 'cellular', 'molecular' ]:
                if topic in self.model:
                    if isinstance(self.model[topic], dict):
                        self.model[topic]['index'] = decoder_index
                        self.model[topic] = self.load_decoder_preset(self.model[topic])
                        self.decoder_by_index.append(self.model[topic])
                        decoder_index +=1

                    elif isinstance(self.model[topic], list):
                        buffer = self.model[topic]
                        self.model[topic] = []
                        for decoder in buffer:
                            decoder['index'] = decoder_index
                            compiled = self.load_decoder_preset(decoder)
                            self.model[topic].append(compiled)
                            self.decoder_by_index.append(compiled)
                            decoder_index +=1

        else: raise BadConfigurationError('unknown preset %s', self.configuration['preset'])

    def load_decoder_preset(self, decoder):
        compiled = decoder
        if 'base' in decoder:
            if 'decoder' in self.configuration and decoder['base'] in self.configuration['decoder']:
                compiled = merge(self.configuration['decoder'][decoder['base']], decoder)
            else:
                raise BadConfigurationError('reference to unknown base decoder {}'.format(decoder['base']))

        if 'unclassified' not in compiled: compiled['unclassified'] = {}
        compiled['unclassified']['index'] = 0
        compiled['unclassified']['description'] = 'unclassified'
        compiled['unclassified']['concentration'] = compiled['noise']
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
        not_noise = 1.0 - compiled['noise']
        for barcode_key in sorted(compiled['codec'].keys()):
            barcode = compiled['codec'][barcode_key]
            barcode['key'] = barcode_key
            barcode['index'] = barcode_index
            barcode['description'] = '{}:{}'.format(barcode_index, '-'.join(barcode['barcode']))
            if 'concentration' not in barcode: barcode['concentration'] = 1.0
            total_concentration += barcode['concentration']

            barcode_length = [ len(i) for i in barcode['barcode'] ]
            if barcode_length != compiled['barcode length']:
                raise BadConfigurationError('found barcode length {} for barcode {} in decoder {} and expected {}'.format(
                    ' '.join([str(i) for i in barcode_length]),
                    barcode_key,
                    str(compiled['index']),
                    ' '.join([str(i) for i in compiled['barcode length']])
                ))

            # convert barcode sequence to a numpy char array
            for segment_index,segment in enumerate(barcode['barcode']):
                barcode['barcode'][segment_index] = numpy.array(list(segment))
            barcode_index += 1

        density = 0
        compiled['barcode distribution.compiled'] = []
        for barcode_key, barcode_value in compiled['codec'].items():
            barcode_value['concentration'] = not_noise * (barcode_value['concentration'] / total_concentration)
            density += barcode_value['concentration']
            compiled['barcode distribution.compiled'].append({'density': density, 'item': barcode_value})

        if compiled['error model reference'] in self.configuration['substitution model']:
            compiled['substitution model'] = { 'frequency': { 'segment': [] } }
            for token_index, token in enumerate(compiled['substitution frequency transform.compiled']['token']):
                compiled['substitution model']['frequency']['segment'].append({'nucleotide substitution frequency': []})

            reference = self.configuration['substitution model'][compiled['error model reference']]
            for token_index, token in enumerate(compiled['substitution frequency transform.compiled']['token']):
                for cycle in range(token['start'], token['end']):
                    compiled['substitution model']['frequency']['segment'][token_index]['nucleotide substitution frequency'].append(deepcopy(reference['segment'][token['segment']]['nucleotide substitution frequency'][cycle]))
        else:
            raise BadConfigurationError('unknown error model reference {}'.format(decoder['error model reference']))

        compiled['quality distribution'] = []
        for q in range(64):
            compiled['quality distribution'].append(0)

        return compiled

    def load_noise_model(self):
        if 'noise model' in self.model:
            if 'reference fasta path' in self.noise_model:
                relative_path = os.path.join(self.instruction['working directory'], self.noise_model['reference fasta path'])
                if os.path.exists(relative_path):
                    with io.open(relative_path, 'rb') as file:
                        content = file.read().decode('ascii').split('\n')
                        if content[0][0] == '>':
                            title = content[0][1:].split(' ')
                            self.noise_model['genome id'] = title[0]
                            self.noise_model['genome comment'] = ' '.join(title[1:])
                            self.noise_model['genome.compiled'] = numpy.array(list(''.join([line.strip() for line in content[1:]])))
                            self.noise_model['nucleotide cardinality'] = len(self.noise_model['genome.compiled'])
                else: raise NoConfigurationFileError('noise model reference fasta {} not found'.format(self.noise_model['reference fasta path']))
            else: raise BadConfigurationError('missing reference fasta path for noise model')

    def load_report(self):
        self.model['count'] = 0
        for decoder in self.decoder_by_index:
            decoder['expected barcode error compensation'] = 0
            decoder['expected barcode error'] = 0
            decoder['unclassified']['count'] = 0

            decoder['simulated nucleotide'] = 0
            decoder['classified count'] = 0

            for barcode in decoder['codec'].values():
                barcode['count'] = 0

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
                    for token, barcode_segment in zip(decoder['transform.compiled']['token'], barcode['barcode']):
                        read['segment'][token['segment']]['fixed'][9][token['start']:token['end']] = barcode_segment

                else:
                    # noise: insert a random piece of PhiX
                    for token, segment_length in zip(decoder['transform.compiled']['token'], decoder['barcode length']):
                        offset = self.sample_random_noise(segment_length)
                        read['segment'][token['segment']]['fixed'][9][token['start']:token['end']] = self.noise_genome[offset:offset + segment_length]

                for token, barcode_segment_index in zip(decoder['transform.compiled']['token'], range(decoder['barcode segment cardinality'])):
                    for token_cycle, segment_cycle in zip(range(token['length']), range(token['start'], token['end'])):
                        nucleotide = read['segment'][token['segment']]['fixed'][9][segment_cycle]
                        quality = decode_phred(read['segment'][token['segment']]['fixed'][10][segment_cycle])
                        decoder['quality distribution'][quality] += 1

                        y = self.phred_scale[quality]['error'] - decoder['expected barcode error compensation']
                        t = decoder['expected barcode error'] + y
                        compensation = (t - decoder['expected barcode error']) - y
                        decoder['expected barcode error'] = t

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
        scale = self.noise_model['nucleotide cardinality'] - size
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
            unclassified['simulated concentration'] = 0
            unclassified['simulated concentration deviation'] = 0
            unclassified['simulated nucleotide'] = unclassified['count'] * decoder['nucleotide cardinality']
            if self.model['count'] > 0:
                unclassified['simulated concentration'] = unclassified['count'] / self.model['count']
            if unclassified['concentration'] > 0:
                unclassified['simulated concentration deviation'] = 1.0 - unclassified['simulated concentration'] / unclassified['concentration']

            # classified
            for barcode in decoder['codec'].values():
                barcode['simulated concentration'] = 0
                barcode['simulated concentration deviation'] = 0

                decoder['classified count'] += barcode['count']
                barcode['simulated nucleotide'] = barcode['count'] * decoder['nucleotide cardinality']
                if self.model['count'] > 0:
                    barcode['simulated concentration'] = barcode['count'] / self.model['count']
                if barcode['concentration'] > 0:
                    barcode['simulated concentration deviation'] = 1.0 - barcode['simulated concentration'] / barcode['concentration']

            del decoder['expected barcode error compensation']
            decoder['simulated noise'] = 0
            decoder['simulated noise deviation'] = 0
            decoder['simulated classified nucleotide'] = decoder['classified count'] * decoder['nucleotide cardinality']
            decoder['simulated nucleotide'] = unclassified['simulated nucleotide'] + decoder['simulated classified nucleotide']
            if decoder['simulated nucleotide'] > 0:
                decoder['expected barcode error'] /= decoder['simulated nucleotide']

            if self.model['count'] > 0:
                decoder['simulated noise'] = unclassified['count'] / self.model['count']
            if decoder['noise'] > 0:
                decoder['simulated noise deviation'] = 1.0 - decoder['simulated noise'] / decoder['noise']
