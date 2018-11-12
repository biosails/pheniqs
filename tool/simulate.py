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

    # Structure of a SAM record
    # -----------------------------------------------------------------------------
    # 0  QNAME   string     Query template NAME
    # 1  FLAG    int        bitwise FLAG
    #    0x1     template having multiple segments in sequencing
    #    0x2     each segment properly aligned according to the aligner
    #    0x4     segment unmapped
    #    0x8     next segment in the template unmapped
    #    0x10    SEQ being reverse complemented
    #    0x20    SEQ of the next segment in the template being reverse complemented
    #    0x40    the first segment in the template
    #    0x80    the last segment in the template
    #    0x100   secondary alignment
    #    0x200   not passing filters, such as platform/vendor quality controls
    #    0x400   PCR or optical duplicate
    #    0x800   supplementary alignment
    # 2  RNAME   string     Reference sequence NAME
    # 3  POS     int        1-based leftmost mapping POSition
    # 4  MAPQ    int        MAPping Quality
    # 5  CIGAR   string     CIGAR string
    # 6  RNEXT   string     Reference name of the mate/next read
    # 7  PNEXT   int        Position of the mate/next read
    # 8  TLEN    int        observed Template LENgth
    # 9  SEQ     string     segment SEQuence
    # 10 QUAL    string     Phred QUALity+33

import numpy.random
from core import *

def print_json(node):
    def handler(o):
        if isinstance(o, numpy.ndarray):
            return ''.join(o)

        elif isinstance(o, Decimal):
            return '{:10.17f}'.format(o)

        return None

    return print(json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4, default=handler))

class Simulate(Pipeline):
    def __init__(self):
        Pipeline.__init__(self, 'simulate')
        self.read = None
        self.input_buffer = None
        self.output_buffer = None
        self.input_handle = None
        self.output_handle = None
        self.pattern = {
            'sam optional field': re.compile('(?P<TAG>[A-Za-z][A-Za-z0-9]):(?P<TYPE>[AifZHB]):(?P<VALUE>[!-~]+)'),
            'sam optional field value': {
                'A': re.compile('^[!-~]$'),
                'i': re.compile('^[-+]?[0-9]+$'),
                'f': re.compile('^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'),
                'Z': re.compile('^[ !-~]+$'),
                'H': re.compile('^[0-9A-F]+$'),
                'B': re.compile('^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+$'),
            }
        }

    @property
    def configuration(self):
        return self.ontology['configuration']

    @property
    def preset(self):
        return self.ontology['preset']

    @property
    def report(self):
        return self.ontology['report']

    @property
    def buffer_capacity(self):
        return self.instruction['capacity']

    @property
    def noise_genome(self):
        return self.preset['noise']['genome']

    def open(self):
        self.input_buffer = []
        self.output_buffer = []

        if self.instruction['input']:
            if os.path.exists(self.instruction['input']):
                self.input_handle = io.open(self.instruction['input'], 'r')
            else:
                raise BadConfigurationError('input {} not found'.format(self.instruction['input']))
        else:
            self.input_handle = sys.stdin

        if self.instruction['output']:
            self.output_handle = io.open(self.instruction['output'], 'w')
        else:
            self.output_handle = sys.stdout

    def close(self):
        if self.input_handle is not None:
            if self.input_handle != sys.stdin:
                self.input_handle.close()
            self.input_handle = None

        if self.output_handle is not None:
            if self.output_handle != sys.stdout:
                self.output_handle.close()
            self.output_handle = None
        Pipeline.close(self)

    def load_configuration(self):
        if os.path.exists(self.instruction['path']):
            self.instruction['working directory'] = os.path.dirname(self.instruction['path'])
            # self.log.debug('loading model %s', self.instruction['path'])
            with io.open(self.instruction['path'], 'rb') as file:
                self.ontology['configuration'] = json.loads(file.read().decode('utf8'))
            self.load_preset()
        else: raise NoConfigurationFileError('configure {} not found'.format(self.instruction['path']))

    def load_preset(self):
        if self.instruction['preset'] in self.configuration['preset']:
            self.ontology['preset'] = deepcopy(self.configuration['preset'][self.instruction['preset']])

            if 'decoder' in self.preset and 'decoder' in self.configuration:
                for decoder_key,decoder in self.preset['decoder'].items():
                    if decoder_key in self.configuration['decoder']:
                        self.preset['decoder'][decoder_key] = merge(self.preset['decoder'][decoder_key], self.configuration['decoder'][decoder_key])

            for decoder_key,decoder in self.preset['decoder'].items():
                decoder['key'] = decoder_key
                decoder['barcode cardinality'] = len(decoder['codec'])
                decoder['transform']['parsed token'] = self.parse_token(decoder['transform']['token'])
                decoder['transform']['parsed error token'] = self.parse_token(decoder['transform']['error token'])

                decoder['barcode length'] = []
                decoder['nucleotide cardinality'] = 0
                for token in decoder['transform']['parsed token']:
                    decoder['barcode length'].append(token['length'])
                    decoder['nucleotide cardinality'] += token['length']
                decoder['barcode segment cardinality'] = len(decoder['barcode length'])

                decoder['error model'] = { 'segment': [] }
                for token_index, token in enumerate(decoder['transform']['parsed error token']):
                    if token['length'] == decoder['transform']['parsed token'][token_index]['length']:
                        decoder['error model']['segment'].append({'error frequency': [], 'substitution frequency': []})
                    else: raise BadConfigurationError('error token {} length does not match decoder'.format(str(token_index)))

                total_concentration = 0
                not_noise = 1.0 - decoder['noise']
                barcode_index = 0
                for barcode_key,barcode_value in decoder['codec'].items():
                    barcode_value['key'] = barcode_key
                    barcode_value['index'] = barcode_index
                    barcode_value['description'] = '{}:{}'.format(barcode_index, '-'.join(barcode_value['barcode']))

                    if 'concentration' not in barcode_value:
                        barcode_value['concentration'] = 1.0
                    total_concentration += barcode_value['concentration']

                    barcode_length = [ len(i) for i in barcode_value['barcode'] ]
                    if barcode_length != decoder['barcode length']:
                        raise BadConfigurationError('found barcode length {} for barcode {} in decoder {} and expected {}'.format(
                            ' '.join([str(i) for i in barcode_length]),
                            barcode_key,
                            decoder_key,
                            ' '.join([str(i) for i in decoder['barcode length']])
                        ))

                    # convert barcode sequence to a numpy char array
                    for segment_index,segment in enumerate(barcode_value['barcode']):
                        barcode_value['barcode'][segment_index] = numpy.array(list(segment))

                    barcode_index += 1

                density = 0
                decoder['distribution'] = []
                for barcode_key,barcode_value in decoder['codec'].items():
                    barcode_value['concentration'] = not_noise * (barcode_value['concentration'] / total_concentration)
                    density += barcode_value['concentration']
                    decoder['distribution'].append({'density': density, 'item': barcode_value})

                if decoder['error model reference'] in self.configuration['error model']:
                    error_model_reference = self.configuration['error model'][decoder['error model reference']]

                    # extract error model tokens
                    for token_index,token in enumerate(decoder['transform']['parsed error token']):
                        for cycle in range(token['start'], token['end']):
                            decoder['error model']['segment'][token_index]['error frequency'].append(deepcopy(error_model_reference['segment'][token['segment']]['error frequency'][cycle]))
                            decoder['error model']['segment'][token_index]['substitution frequency'].append(deepcopy(error_model_reference['segment'][token['segment']]['substitution frequency'][cycle]))

                    # normalize error frequencies on the tokenized error model
                    decoder['error model']['barcode expected error'] = 0
                    for segment in decoder['error model']['segment']:
                        for cycle in segment['error frequency']:
                            cycle_expected_error = sum(cycle.values()) / 4.0
                            decoder['error model']['barcode expected error'] += cycle_expected_error

                    decoder['error model']['barcode expected error per nucleotide'] = decoder['error model']['barcode expected error'] / decoder['nucleotide cardinality']
                    decoder['error model']['barcode expected error factor'] = decoder['expected error frequency'] / decoder['error model']['barcode expected error per nucleotide']

                    for segment in decoder['error model']['segment']:
                        segment['nomrlaized error frequency'] = []
                        for cycle in segment['error frequency']:
                            normalized_cycle = {}
                            for k,v in cycle.items():
                                normalized_cycle[k] = v * decoder['error model']['barcode expected error factor']
                            segment['nomrlaized error frequency'].append(normalized_cycle)

                    decoder['error model']['nomrlaized barcode expected error'] = 0
                    for segment in decoder['error model']['segment']:
                        for cycle in segment['nomrlaized error frequency']:
                            cycle_expected_error = sum(cycle.values()) / 4.0
                            decoder['error model']['nomrlaized barcode expected error'] += cycle_expected_error
                    decoder['error model']['nomrlaized barcode expected error per nucleotide'] = decoder['error model']['nomrlaized barcode expected error'] / decoder['nucleotide cardinality']

                else: raise BadConfigurationError('unknown error model reference {}'.format(decoder['error model reference']))

            self.load_noise_model()
            self.load_report()

        else: raise BadConfigurationError('unknown preset %s', self.configuration['preset'])

    def load_noise_model(self):
        if 'noise' in self.preset:
            if 'reference fasta path' in self.preset['noise']:
                relative_path = os.path.join(self.instruction['working directory'], self.preset['noise']['reference fasta path'])
                if os.path.exists(relative_path):
                    with io.open(relative_path, 'rb') as file:
                        content = file.read().decode('ascii').split('\n')
                        if content[0][0] == '>':
                            title = content[0][1:].split(' ')
                            self.preset['noise']['id'] = title[0]
                            self.preset['noise']['comment'] = ' '.join(title[1:])
                            self.preset['noise']['genome'] = numpy.array(list(''.join([line.strip() for line in content[1:]])))
                            self.preset['noise']['nucleotide cardinality'] = len(self.preset['noise']['genome'])
                else: raise NoConfigurationFileError('noise model reference fasta {} not found'.format(self.preset['noise']['reference fasta path']))
            else: raise BadConfigurationError('missing reference fasta path for noise model')

    def load_report(self):
        self.ontology['report'] = {
            'count': 0,
            'decoder': {}
        }
        for decoder_key, decoder in self.preset['decoder'].items():
            self.report['decoder'][decoder_key] = {
                'simulated nucleotide': 0,
                'simulated noise substitution': 0,
                'simulated classified substitution': 0,
                'simulated substitution': 0,
                'classified count': 0,
                'noise count': 0,
                'barcode cardinality': deepcopy(decoder['barcode cardinality']),
                'barcode length': deepcopy(decoder['barcode length']),
                'expected error frequency': deepcopy(decoder['expected error frequency']),
                'noise': deepcopy(decoder['noise']),
                'nucleotide cardinality': deepcopy(decoder['nucleotide cardinality']),
                'codec': deepcopy(decoder['codec']),
                'error model': {
                    'barcode expected error': deepcopy(decoder['error model']['nomrlaized barcode expected error']),
                    'barcode expected error per nucleotide': deepcopy(decoder['error model']['nomrlaized barcode expected error per nucleotide']),
                    'segment': []
                },
            }
            for barcode in self.report['decoder'][decoder_key]['codec'].values():
                del barcode['key']
                del barcode['description']
                barcode['count'] = 0
                barcode['simulated substitution'] = 0

            for segment in decoder['error model']['segment']:
                self.report['decoder'][decoder_key]['error model']['segment'].append({
                    'error frequency': deepcopy(segment['nomrlaized error frequency']),
                    'substitution frequency': deepcopy(segment['substitution frequency'])
                })

    def execute(self):
        if self.action == 'barcode':
            self.execute_barcode()

    def execute_barcode(self):
        self.load_configuration()

        self.open()
        while(self.replenish()):
            self.manipulate()
            self.flush()
        self.close()

        self.finalize()
        print_json(self.report)

    def replenish(self):
        for line in self.input_handle:
            if line:
                if line[0] =='@':
                    # redirect header line
                    self.output_handle.write(line)
                else:
                    segment = {
                        'fixed': line.strip().split('\t'),
                        'auxiliary': {}
                    }
                    if len(segment['fixed']) > 10:
                        # parse AUX
                        if len(segment['fixed']) > 11:
                            for item in segment['fixed'][11:]:
                                field = item.split(':')
                                if len(field) == 3:
                                    tag = { 'TAG': field[0], 'TYPE': field[1], 'VALUE': field[2] }
                                    # if self.pattern['sam optional field value'][tag['TYPE']].match(tag['VALUE']):
                                    if tag['TYPE'] is 'i':
                                        tag['VALUE'] = int(tag['VALUE'])
                                    elif tag['TYPE'] is 'f':
                                        tag['VALUE'] = float(tag['VALUE'])
                                    elif tag['TYPE'] is 'H':
                                        raise UnsupportedError('hex format byte array not supported')
                                    elif tag['TYPE'] is 'B':
                                        raise UnsupportedError('numeric array not supported')
                                    segment['auxiliary'][tag['TAG']] = tag
                                    # else: log.error('ignoring invalid %s optional field %s', tag['TYPE'], tag['VALUE'])
                                else:
                                    log.error('ignoring invalid auxiliary tag %s', field)

                        segment['fixed'][9] = numpy.array(list(segment['fixed'][9]))

                        if self.read:
                            # check QNAME is the same as the buffer
                            if self.read[0]['fixed'][0] == segment['fixed'][0]:
                                self.read.append(segment)
                            else:
                                # new read, write the complete one to the buffer and start over
                                self.input_buffer.append(self.read)
                                self.read = [ segment ]
                                if len(self.input_buffer) >= self.buffer_capacity: break
                        else:
                            self.read = [ segment ]
                    else:
                        self.log.error('invalid sam syntax %s', line)
                        raise SequenceError(line)

        if len(self.input_buffer) == 0 and self.read:
            self.input_buffer.append(self.read)
            self.read = None

        self.report['count'] += len(self.input_buffer)
        return len(self.input_buffer)

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)

            for decoder in self.preset['decoder'].values():
                barcode = self.pick_barcode(decoder)
                if barcode:
                    self.report['decoder'][decoder['key']]['codec'][barcode['key']]['count'] += 1
                    for segment in read:
                        segment['auxiliary'][decoder['true barcode tag']] = {
                            'TAG': decoder['true barcode tag'],
                            'TYPE': 'Z',
                            'VALUE': '-'.join([''.join(s) for s in barcode['barcode']])
                        }

                    for token,barcode_segment_index,barcode_segment in zip(decoder['transform']['parsed token'], range(decoder['barcode segment cardinality']), barcode['barcode']):
                        read[token['segment']]['fixed'][9][token['start']:token['end']] = barcode_segment
                        for token_cycle,segment_cycle in zip(range(token['length']), range(token['start'], token['end'])):
                            nucleotide = read[token['segment']]['fixed'][9][segment_cycle]
                            if self.error_occured(decoder, barcode_segment_index, token_cycle, nucleotide):
                                self.report['decoder'][decoder['key']]['codec'][barcode['key']]['simulated substitution'] += 1
                                substitute = self.substitute_with(decoder, barcode_segment_index, token_cycle, nucleotide)
                                read[token['segment']]['fixed'][9][segment_cycle] = substitute
                                # self.log.debug('substitution event %s : %s -> %s', cycle, nucleotide, substitute)
                else:
                    # noise: insert a random piece of PhiX
                    self.report['decoder'][decoder['key']]['noise count'] += 1
                    for token,barcode_segment_index,barcode_segment_length in zip(decoder['transform']['parsed token'], range(decoder['barcode segment cardinality']), decoder['barcode length']):
                        offset = self.sample_random_noise(barcode_segment_length)
                        read[token['segment']]['fixed'][9][token['start']:token['end']] = self.noise_genome[offset:offset + barcode_segment_length]
                        for token_cycle,segment_cycle in zip(range(token['length']), range(token['start'], token['end'])):
                            nucleotide = read[token['segment']]['fixed'][9][segment_cycle]
                            if self.error_occured(decoder, barcode_segment_index, token_cycle, nucleotide):
                                self.report['decoder'][decoder['key']]['simulated noise substitution'] += 1
                                substitute = self.substitute_with(decoder, barcode_segment_index, token_cycle, nucleotide)
                                read[token['segment']]['fixed'][9][segment_cycle] = substitute

            self.output_buffer.append(read)

    def flush(self):
        while(len(self.output_buffer) > 0):
            read = self.output_buffer.pop(0)
            for segment in read:
                fixed = '\t'.join([
                    segment['fixed'][0],
                    segment['fixed'][1],
                    segment['fixed'][2],
                    segment['fixed'][3],
                    segment['fixed'][4],
                    segment['fixed'][5],
                    segment['fixed'][6],
                    segment['fixed'][7],
                    segment['fixed'][8],
                    ''.join(segment['fixed'][9]),
                    segment['fixed'][10]
                ])
                if segment['auxiliary']:
                    auxiliary = []
                    for tag in segment['auxiliary'].values():
                        value = None
                        if tag['TYPE'] is 'A':
                            value = tag['VALUE']
                        elif tag['TYPE'] is 'i':
                            value = str(tag['VALUE'])
                        elif tag['TYPE'] is 'f':
                            value = str(tag['VALUE'])
                        elif tag['TYPE'] is 'Z':
                            value = tag['VALUE']
                        elif tag['TYPE'] is 'H':
                            raise UnsupportedError('hex format byte array not supported')
                        elif tag['TYPE'] is 'B':
                            raise UnsupportedError('numeric array not supported')
                        auxiliary.append('{}:{}:{}'.format(tag['TAG'], tag['TYPE'], value))
                    auxiliary = '\t'.join(auxiliary)
                self.output_handle.write('{}\t{}\n'.format(fixed, auxiliary))

    def finalize(self):
        for decoder in self.report['decoder'].values():
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
                if self.report['count'] > 0:
                    barcode['simulated concentration'] = barcode['count'] / self.report['count']
                if barcode['concentration'] > 0:
                    barcode['simulated concentration deviation'] = 1.0 - barcode['simulated concentration'] / barcode['concentration']

            decoder['substitution rate per nucleotide'] = 0
            decoder['substitution rate per classified nucleotide'] = 0
            decoder['substitution rate per noise nucleotide'] = 0
            decoder['simulated noise'] = 0
            decoder['simulated noise deviation'] = 0

            decoder['simulated classified nucleotide'] = decoder['classified count'] * decoder['nucleotide cardinality']
            decoder['simulated noise nucleotide'] = decoder['noise count'] * decoder['nucleotide cardinality']
            decoder['simulated nucleotide'] = decoder['simulated classified nucleotide'] + decoder['simulated noise nucleotide']
            decoder['simulated substitution'] += decoder['simulated classified substitution'] + decoder['simulated noise substitution']
            if decoder['simulated nucleotide'] > 0:
                decoder['substitution rate per nucleotide'] = decoder['simulated substitution'] / decoder['simulated nucleotide']
            if decoder['simulated classified nucleotide'] > 0:
                decoder['substitution rate per classified nucleotide'] = decoder['simulated classified substitution'] / decoder['simulated classified nucleotide']
            if decoder['simulated noise nucleotide'] > 0:
                decoder['substitution rate per noise nucleotide'] = decoder['simulated noise substitution'] / decoder['simulated noise nucleotide']
            if self.report['count'] > 0:
                decoder['simulated noise'] = decoder['noise count'] / self.report['count']
            if decoder['noise'] > 0:
                decoder['simulated noise deviation'] = 1.0 - decoder['simulated noise'] / decoder['noise']

    def parse_token(self, pattern_array):
        parsed_pattern_array = []
        for pattern in pattern_array:
            parsed = pattern.split(':')
            token = { 'segment': parsed[0], 'start': parsed[1], 'end': parsed[2] }

            if len(token['segment']) > 0:
                token['segment'] = int(token['segment'])
            else:
                raise BadConfigurationError('segment index is missing in {}'.format(pattern))

            if len(token['start']) is 0:
                token['start'] = 0
            token['start'] = int(token['start'])

            if len(token['end']) > 0:
                token['end'] = int(token['end'])
            else:
                raise BadConfigurationError('end position is missing in {}'.format(pattern))

            token['length'] = token['end'] - token['start']

            parsed_pattern_array.append(token)

        return parsed_pattern_array

    def sample_random_noise(self, size):
        scale = self.preset['noise']['nucleotide cardinality'] - size
        event = numpy.random.random()
        return int(event * scale)

    def pick_barcode(self, decoder):
        barcode = None
        event = numpy.random.random()
        for option in decoder['distribution']:
            if event < option['density']:
                barcode = option['item']
                # self.log.debug('barcode %s picked by %.15f', barcode['description'], event)
                break
        # if barcode is None: self.log.debug('noise picked %.15f', event)
        return barcode

    def error_occured(self, decoder, segment, cycle, nucleotide):
        if nucleotide != 'N':
            event = numpy.random.random()
            return event < decoder['error model']['segment'][segment]['nomrlaized error frequency'][cycle][nucleotide]
        else:
            return False

    def substitute_with(self, decoder, segment, cycle, nucleotide):
        event = numpy.random.random()
        for option in decoder['error model']['segment'][segment]['substitution frequency'][cycle][nucleotide]:
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
