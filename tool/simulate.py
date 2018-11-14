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

from sam import *

class Simulate(SamTranscode):
    def __init__(self):
        SamTranscode.__init__(self, 'simulate')
        self.ontology['barcode distribution'] = {}

    @property
    def configuration(self):
        return self.ontology['configuration']

    @property
    def barcode_distribution(self):
        return self.ontology['barcode distribution']

    @property
    def noise_model(self):
        return self.preset['noise model']

    @property
    def noise_genome(self):
        return self.preset['noise model']['genome']

    @property
    def read_count(self):
        return self.ontology['preset']['count']

    def load(self):
        self.load_configuration()
        self.load_preset()

    def load_configuration(self):
        if os.path.exists(self.instruction['path']):
            self.instruction['working directory'] = os.path.dirname(self.instruction['path'])
            self.log.debug('loading model %s', self.instruction['path'])
            with io.open(self.instruction['path'], 'rb') as file:
                self.ontology['configuration'] = json.loads(file.read().decode('utf8'))
        else: raise NoConfigurationFileError('configure {} not found'.format(self.instruction['path']))

    def load_preset(self):
        if self.instruction['preset'] in self.configuration['preset']:
            self.ontology['preset'] = deepcopy(self.configuration['preset'][self.instruction['preset']])
            self.ontology['preset']['run count'] = str(self.instruction['number'])

            if 'decoder' in self.preset and 'decoder' in self.configuration:
                for decoder_key,decoder in self.preset['decoder'].items():
                    decoder['unclassified'] = {
                        'index': 0,
                        'concentration': decoder['noise'],
                        'description': 'unclassified',
                    }
                    if decoder_key in self.configuration['decoder']:
                        self.preset['decoder'][decoder_key] = merge(self.preset['decoder'][decoder_key], self.configuration['decoder'][decoder_key])
                        self.preset['decoder'][decoder_key]['expected substitution frequency'] = self.instruction['error']
            decoder_index = 0
            for key,decoder in self.preset['decoder'].items():
                decoder['index'] = decoder_index
                decoder['key'] = key
                decoder['barcode length'] = []
                decoder['nucleotide cardinality'] = 0
                decoder['barcode cardinality'] = len(decoder['codec'])

                self.token_model[key] = {
                    'token': self.parse_token(decoder['transform']['token']),
                    'error token': self.parse_token(decoder['transform']['error token']),
                }
                for token in self.token_model[key]['token']:
                    decoder['barcode length'].append(token['length'])
                    decoder['nucleotide cardinality'] += token['length']

                decoder['barcode segment cardinality'] = len(decoder['barcode length'])
                if len(self.token_model[key]['error token']) != decoder['barcode segment cardinality']:
                    raise BadConfigurationError('error token array size does not match decoder {} token array'.format(key))

                # assign index to each barcode and normalize the prior
                # barcode index 0 is reserved for unclassified
                barcode_index = 1
                total_concentration = 0
                not_noise = 1.0 - decoder['noise']
                for barcode_key in sorted(decoder['codec'].keys()):
                    barcode = decoder['codec'][barcode_key]
                    barcode['key'] = barcode_key
                    barcode['index'] = barcode_index
                    barcode['description'] = '{}:{}'.format(barcode_index, '-'.join(barcode['barcode']))
                    if 'concentration' not in barcode: barcode['concentration'] = 1.0
                    total_concentration += barcode['concentration']

                    barcode_length = [ len(i) for i in barcode['barcode'] ]
                    if barcode_length != decoder['barcode length']:
                        raise BadConfigurationError('found barcode length {} for barcode {} in decoder {} and expected {}'.format(
                            ' '.join([str(i) for i in barcode_length]),
                            barcode_key,
                            decoder_key,
                            ' '.join([str(i) for i in decoder['barcode length']])
                        ))

                    # convert barcode sequence to a numpy char array
                    for segment_index,segment in enumerate(barcode['barcode']):
                        barcode['barcode'][segment_index] = numpy.array(list(segment))
                    barcode_index += 1

                density = 0
                self.barcode_distribution[key] = []
                for barcode_key,barcode_value in decoder['codec'].items():
                    barcode_value['concentration'] = not_noise * (barcode_value['concentration'] / total_concentration)
                    density += barcode_value['concentration']
                    self.barcode_distribution[key].append({'density': density, 'item': barcode_value})

                # load a substitution model
                if decoder['error model reference'] in self.configuration['substitution model']:
                    model = { 'segment': [] }
                    for token_index, token in enumerate(self.token_model[key]['error token']):
                        if token['length'] == self.token_model[key]['token'][token_index]['length']:
                            model['segment'].append({'substitution occurrence frequency': [], 'nucleotide substitution frequency': []})
                        else: raise BadConfigurationError('error token {} length does not match decoder'.format(str(token_index)))

                    reference = self.configuration['substitution model'][decoder['error model reference']]
                    for token_index,token in enumerate(self.token_model[key]['error token']):
                        for cycle in range(token['start'], token['end']):
                            model['segment'][token_index]['substitution occurrence frequency'].append(deepcopy(reference['segment'][token['segment']]['substitution occurrence frequency'][cycle]))
                            model['segment'][token_index]['nucleotide substitution frequency'].append(deepcopy(reference['segment'][token['segment']]['nucleotide substitution frequency'][cycle]))

                    # derive a normalized substitution model
                    normalized = {
                        'segment': [],
                        'barcode expected error': 0,
                        'barcode expected error per nucleotide': 0,
                    }
                    for token_index, token in enumerate(self.token_model[key]['error token']):
                        normalized['segment'].append({'substitution occurrence frequency': [], 'nucleotide substitution frequency': []})

                    barcode_expected_error = 0
                    for segment in model['segment']:
                        for cycle in segment['substitution occurrence frequency']:
                            barcode_expected_error += (sum(cycle.values()) / 4.0)
                    barcode_expected_error_per_nucleotide = barcode_expected_error / decoder['nucleotide cardinality']

                    substitution_factor = decoder['expected substitution frequency'] / barcode_expected_error_per_nucleotide
                    for model_segment, normalized_segment in zip(model['segment'], normalized['segment']):
                        normalized_segment['nucleotide substitution frequency'] = deepcopy(model_segment['nucleotide substitution frequency'])
                        for model_cycle in model_segment['substitution occurrence frequency']:
                            normalized_cycle = {}
                            for k,v in model_cycle.items():
                                normalized_cycle[k] = v * substitution_factor
                            normalized_segment['substitution occurrence frequency'].append(normalized_cycle)

                    for segment in normalized['segment']:
                        for cycle in segment['substitution occurrence frequency']:
                            normalized['barcode expected error'] += (sum(cycle.values()) / 4.0)
                    normalized['barcode expected error per nucleotide'] = normalized['barcode expected error'] / decoder['nucleotide cardinality']
                    decoder['substitution model'] = normalized

                    self.log.info('loaded %s decoder', key)
                else: raise BadConfigurationError('unknown error model reference {}'.format(decoder['error model reference']))

                self.decoder_by_index.append(decoder)
                decoder_index += 1

            self.load_noise_model()
            self.load_report()

        else: raise BadConfigurationError('unknown preset %s', self.configuration['preset'])

    def load_noise_model(self):
        if 'reference fasta path' in self.noise_model:
            relative_path = os.path.join(self.instruction['working directory'], self.noise_model['reference fasta path'])
            if os.path.exists(relative_path):
                with io.open(relative_path, 'rb') as file:
                    content = file.read().decode('ascii').split('\n')
                    if content[0][0] == '>':
                        title = content[0][1:].split(' ')
                        self.noise_model['genome id'] = title[0]
                        self.noise_model['genome comment'] = ' '.join(title[1:])
                        self.noise_model['genome'] = numpy.array(list(''.join([line.strip() for line in content[1:]])))
                        self.noise_model['nucleotide cardinality'] = len(self.noise_model['genome'])
            else: raise NoConfigurationFileError('noise model reference fasta {} not found'.format(self.noise_model['reference fasta path']))
        else: raise BadConfigurationError('missing reference fasta path for noise model')

    def load_report(self):
        self.preset['count'] = 0
        for decoder in self.decoder_by_index:
            decoder['unclassified']['count'] = 0
            decoder['unclassified']['simulated substitution'] = 0

            decoder['simulated nucleotide'] = 0
            decoder['simulated classified substitution'] = 0
            decoder['simulated substitution'] = 0
            decoder['classified count'] = 0

            for barcode in decoder['codec'].values():
                barcode['count'] = 0
                barcode['simulated substitution'] = 0

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            self.preset['count'] += 1

            # convert segment SEQ to a numpy char array
            for segment in read['segment']:
                segment['fixed'][9] = numpy.array(list(segment['fixed'][9]))

            read['index'] = self.read_count
            read['fact'] = { 'sample': [], 'cellular': [] }
            for decoder in self.decoder_by_index:

                barcode = self.pick_barcode(decoder)
                barcode['count'] += 1

                if decoder['type'] == 'sample':
                    read['fact']['sample'].append(str(barcode['index']))
                elif decoder['type'] == 'cellular':
                    read['fact']['cellular'].append(str(barcode['index']))

                if barcode['index'] > 0:
                    for token,barcode_segment_index,barcode_segment in zip(self.token_model[decoder['key']]['token'], range(decoder['barcode segment cardinality']), barcode['barcode']):
                        read['segment'][token['segment']]['fixed'][9][token['start']:token['end']] = barcode_segment
                        for token_cycle,segment_cycle in zip(range(token['length']), range(token['start'], token['end'])):
                            nucleotide = read['segment'][token['segment']]['fixed'][9][segment_cycle]
                            if self.error_occured(decoder, barcode_segment_index, token_cycle, nucleotide):
                                barcode['simulated substitution'] += 1
                                substitute = self.substitute_with(decoder, barcode_segment_index, token_cycle, nucleotide)
                                read['segment'][token['segment']]['fixed'][9][segment_cycle] = substitute
                else:
                    # noise: insert a random piece of PhiX
                    for token,barcode_segment_index,barcode_segment_length in zip(self.token_model[decoder['key']]['token'], range(decoder['barcode segment cardinality']), decoder['barcode length']):
                        offset = self.sample_random_noise(barcode_segment_length)
                        read['segment'][token['segment']]['fixed'][9][token['start']:token['end']] = self.noise_genome[offset:offset + barcode_segment_length]
                        for token_cycle,segment_cycle in zip(range(token['length']), range(token['start'], token['end'])):
                            nucleotide = read['segment'][token['segment']]['fixed'][9][segment_cycle]
                            if self.error_occured(decoder, barcode_segment_index, token_cycle, nucleotide):
                                barcode['simulated substitution'] += 1
                                substitute = self.substitute_with(decoder, barcode_segment_index, token_cycle, nucleotide)
                                read['segment'][token['segment']]['fixed'][9][segment_cycle] = substitute

            self.simulate_qname(read)

            # convert segment SEQ back to string
            for segment in read['segment']:
                segment['fixed'][9] = ''.join(segment['fixed'][9])

            self.output_buffer.append(read)

    def finalize(self):
        for decoder in self.preset['decoder'].values():
            # unclassified
            unclassified = decoder['unclassified']
            unclassified['substitution rate per barcode'] = 0
            unclassified['substitution rate per nucleotide'] = 0
            unclassified['simulated concentration'] = 0
            unclassified['simulated concentration deviation'] = 0
            if unclassified['count'] > 0:
                unclassified['substitution rate per barcode'] = unclassified['simulated substitution'] / unclassified['count']
            unclassified['simulated nucleotide'] = unclassified['count'] * decoder['nucleotide cardinality']
            if unclassified['simulated nucleotide'] > 0:
                unclassified['substitution rate per nucleotide'] = unclassified['simulated substitution'] / unclassified['simulated nucleotide']
            if self.preset['count'] > 0:
                unclassified['simulated concentration'] = unclassified['count'] / self.preset['count']
            if unclassified['concentration'] > 0:
                unclassified['simulated concentration deviation'] = 1.0 - unclassified['simulated concentration'] / unclassified['concentration']

            # classified
            for barcode in decoder['codec'].values():
                barcode['substitution rate per barcode'] = 0
                barcode['substitution rate per nucleotide'] = 0
                barcode['simulated concentration'] = 0
                barcode['simulated concentration deviation'] = 0

                decoder['classified count'] += barcode['count']
                decoder['simulated classified substitution'] += barcode['simulated substitution']
                if barcode['count'] > 0:
                    barcode['substitution rate per barcode'] = barcode['simulated substitution'] / barcode['count']
                barcode['simulated nucleotide'] = barcode['count'] * decoder['nucleotide cardinality']
                if barcode['simulated nucleotide'] > 0:
                    barcode['substitution rate per nucleotide'] = barcode['simulated substitution'] / barcode['simulated nucleotide']
                if self.preset['count'] > 0:
                    barcode['simulated concentration'] = barcode['count'] / self.preset['count']
                if barcode['concentration'] > 0:
                    barcode['simulated concentration deviation'] = 1.0 - barcode['simulated concentration'] / barcode['concentration']

            decoder['substitution rate per nucleotide'] = 0
            decoder['substitution rate per classified nucleotide'] = 0
            decoder['unclassified']['substitution rate per nucleotide'] = 0
            decoder['simulated noise'] = 0
            decoder['simulated noise deviation'] = 0

            decoder['simulated classified nucleotide'] = decoder['classified count'] * decoder['nucleotide cardinality']
            decoder['simulated nucleotide'] = unclassified['simulated nucleotide'] + decoder['simulated classified nucleotide']
            decoder['simulated substitution'] += unclassified['simulated substitution'] + decoder['simulated classified substitution']
            if decoder['simulated nucleotide'] > 0:
                decoder['substitution rate per nucleotide'] = decoder['simulated substitution'] / decoder['simulated nucleotide']
            if decoder['simulated classified nucleotide'] > 0:
                decoder['substitution rate per classified nucleotide'] = decoder['simulated classified substitution'] / decoder['simulated classified nucleotide']

            if self.preset['count'] > 0:
                decoder['simulated noise'] = unclassified['count'] / self.preset['count']
            if decoder['noise'] > 0:
                decoder['simulated noise deviation'] = 1.0 - decoder['simulated noise'] / decoder['noise']

    def simulate_qname(self, read):
        # @HWI-ST911:232:HABDFADXX:1:1101:1224:1932 1:N:0:CGATGT
        # name    0:1:2:3:4:5:6
        # 0   Instrument ID   HWI-ST911     PHENIQS
        # 1   Run count       232           1
        # 2   Flowcell ID     HABDFADXX     simulated flowcell id
        # 3   Lane number     1             is noise
        # 4   Tile number     1101          Sample barcode index
        # 5   x coordinate    1224          First Cellular barcode index
        # 6   y coordinate    1932          Second Cellular barcode index
        #
        # comment 7:8:9:10
        # 7   Segment number  1               uint8_t
        # 8   Filtered        N               N|Y
        # 9   control number  0               uint16_t
        # 10  Barcode         CGATGT          char*
        qname = [ self.instrument_id, self.run_count, self.flowcell_id ]
        qname.append(str(read['index']).rjust(10,'0'))

        if 'sample' in read['fact']:
            qname.append('S{}'.format('-'.join(read['fact']['sample'])))
        else:
            qname.append('S')

        if 'cellular' in read['fact']:
            qname.append('C{}'.format('-'.join(read['fact']['cellular'])))
        else:
            qname.append('C')

        simulated_qname = ':'.join(qname)
        for segment in read['segment']:
            segment['fixed'][0] = simulated_qname

    def sample_random_noise(self, size):
        scale = self.preset['noise model']['nucleotide cardinality'] - size
        event = numpy.random.random()
        return int(event * scale)

    def pick_barcode(self, decoder):
        barcode = None
        event = numpy.random.random()
        for option in self.barcode_distribution[decoder['key']]:
            if event < option['density']:
                barcode = option['item']
                self.log.debug('barcode %s picked by %.15f', barcode['description'], event)
                break
        if barcode is None:
            barcode = decoder['unclassified']
            self.log.debug('unclassified picked %.15f', event)
        return barcode

    def error_occured(self, decoder, segment, cycle, nucleotide):
        if nucleotide != 'N':
            event = numpy.random.random()
            return event < decoder['substitution model']['segment'][segment]['substitution occurrence frequency'][cycle][nucleotide]
        else:
            return False

    def substitute_with(self, decoder, segment, cycle, nucleotide):
        event = numpy.random.random()
        for option in decoder['substitution model']['segment'][segment]['nucleotide substitution frequency'][cycle][nucleotide]:
            if event < option['density']:
                return option['to']
        return 'N'

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)

    pipeline = None

    try:
        pipeline = Simulate()
        pipeline.execute()

    except (
        PermissionDeniedError,
        NoOverwriteError,
        DownloadError,
        CommandFailedError,
        NoConfigurationFileError,
        BadConfigurationError,
        UnsupportedError,
        SequenceError
    ) as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except(KeyboardInterrupt, SystemExit) as e:
        if e.code != 0:
            logging.getLogger('main').critical(e)
            sys.exit(1)

    finally:
        if pipeline: pipeline.close()

    sys.exit(0)

if __name__ == '__main__':
    main()
