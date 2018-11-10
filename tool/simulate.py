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

def to_json(node):
    def handler(o):
        result = None
        if isinstance(o, numpy.ndarray):
            result = ''.join(o)
        return result
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4, default=handler)

class Simulate(Pipeline):
    def __init__(self):
        Pipeline.__init__(self, 'simulate')
        self.model = None
        self.preset = None
        self.read = None
        self.input_handle = sys.stdin
        self.output_handle = sys.stdout
        self.input_buffer = []
        self.output_buffer = []
        self.pattern = {
            'sam record': re.compile(
                r"""^
                (?P<QNAME>[!-?A-~]{1,254})\t
                (?P<FLAG>[0-9]+)\t
                (?P<RNAME>\*|[!-()+-<>-~][!-~]*)\t
                (?P<POS>[0-9]+)\t
                (?P<MAPQ>[0-9]+)\t
                (?P<CIGAR>\*|([0-9]+[MIDNSHPX=])+)\t
                (?P<RNEXT>\*|=|[!-()+-<>-~][!-~]*)\t
                (?P<PNEXT>[0-9]+)\t
                (?P<TLEN>[-+]?[0-9]+)\t
                (?P<SEQ>\*|[A-Za-z=.]+)\t
                (?P<QUAL>[!-~]*)
                (?P<AUX>(\t([A-Za-z][A-Za-z0-9]:[AifZHB]:[^\s]+))+)?
                $""",
                re.VERBOSE
            ),
            'sam optional field': re.compile('(?P<TAG>[A-Za-z][A-Za-z0-9]):(?P<TYPE>[AifZHB]):(?P<VALUE>[!-~]+)'),
            'sam optional field value': {
                'A': re.compile('^[!-~]$'),
                'i': re.compile('^[-+]?[0-9]+$'),
                'f': re.compile('^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'),
                'Z': re.compile('^[ !-~]+$'),
                'H': re.compile('^[0-9A-F]+$'),
                'B': re.compile('^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+$'),
            },
            'sam record field': [
                'QNAME',
                'FLAG',
                'RNAME',
                'POS',
                'MAPQ',
                'CIGAR',
                'RNEXT',
                'PNEXT',
                'TLEN',
                'SEQ',
                'QUAL',
                'AUX',
            ],
            'sam auxiliary encode': '{TAG}:{TYPE}:{VALUE}',
        }
        self.pattern['sam record encode']: '\t'.join(['{{{}}}'.format(i) for i in self.pattern['sam record field']])

    @property
    def buffer_capacity(self):
        return self.instruction['capacity']

    @property
    def error_frequency(self):
        return self.instruction['frequency']

    def execute(self):
        if self.action == 'barcode':
            self.execute_barcode()

    def execute_barcode(self):
        self.load_model()
        self.load_preset()
        while(self.replenish()):
            self.manipulate()
            self.flush()

    def load_preset(self):
        self.preset = deepcopy(self.model['preset'][self.instruction['preset']])
        for k,v in self.preset.items():
            if k in self.model['vocabulary']:
                self.preset[k] = merge(self.preset[k], self.model['vocabulary'][k])

        for k,v in self.preset.items():
            v['key'] = k
            v['barcode cardinality'] = len(v['codec'])
            v['barcode length'] = None
            v['transform']['parsed token'] = []
            for pattern in v['transform']['token']:
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

                v['transform']['parsed token'].append(token)

            barcode_index = 0
            total_prior = 0
            not_noise = 1.0 - v['noise']
            for barcode_key, barcode_value in v['codec'].items():
                barcode_value['index'] = barcode_index
                if 'concentration' not in barcode_value:
                    barcode_value['concentration'] = 1
                total_prior += barcode_value['concentration']
                if v['barcode length'] is None:
                    v['barcode length'] = [ len(i) for i in barcode_value['barcode'] ]
                else:
                    l = [ len(i) for i in barcode_value['barcode'] ]
                    if l != v['barcode length']:
                        raise BadConfigurationError('barcode {} has wrong length {}'.format(barcode_key, ' '.join([str(i) for i in l])))
                barcode_index += 1

                for segment_index,segment in enumerate(barcode_value['barcode']):
                    barcode_value['barcode'][segment_index] = numpy.array(list(segment))

            density = 0
            v['distribution'] = []
            for barcode_key, barcode_value in v['codec'].items():
                barcode_value['concentration'] = not_noise * (barcode_value['concentration'] / total_prior)
                density += barcode_value['concentration']
                v['distribution'].append({'density': density, 'item': barcode_value})
        # print(to_json(self.preset))

    def load_model(self):
        if os.path.exists(self.instruction['model']):
            self.log.debug('loading model %s', self.instruction['model'])
            with io.open(self.instruction['model'], 'rb') as file:
                self.model = json.loads(file.read().decode('utf8'))

            error_rate = 0
            for position in self.model['error frequency']:
                error_rate += (sum(position.values()) / len(position))
            error_rate /= len(self.model['error frequency'])
            self.model['error frequency per base'] = error_rate
            self.model['error occurance factor'] = self.error_frequency / self.model['error frequency per base']
            # print(to_json(self.model))

        else:
            raise NoConfigurationFileError('configure {} not found'.format(self.instruction['model']))

    def replenish(self):
        for line in self.input_handle:
            if line:
                if line[0] =='@':
                    pass
                    # header line
                    # self.output_handle.write(line)
                else:
                    match = self.pattern['sam record'].search(line)
                    if match:
                        segment = match.groupdict()
                        for field in [
                            'FLAG',
                            'POS',
                            'MAPQ',
                            'PNEXT',
                            'TLEN',
                        ]:
                            segment[field] = int(segment[field])

                        segment['SEQ'] = numpy.array(list(segment['SEQ']))
                        # segment['QUAL'] = numpy.array(list(segment['QUAL']))
                        if 'AUX' in segment:
                            optional = segment['AUX'].strip('\t').split('\t')
                            auxiliary = {}
                            for o in optional:
                                m = self.pattern['sam optional field'].search(o)
                                if m:
                                    tag = m.groupdict()
                                    if self.pattern['sam optional field value'][tag['TYPE']].match(tag['VALUE']):
                                        if tag['TYPE'] is 'i':
                                            tag['VALUE'] = int(tag['VALUE'])
                                        elif tag['TYPE'] is 'f':
                                            tag['VALUE'] = float(tag['VALUE'])
                                        elif tag['TYPE'] is 'H':
                                            raise UnsupportedError('hex format byte array not supported')
                                        elif tag['TYPE'] is 'B':
                                            raise UnsupportedError('numeric array not supported')
                                        auxiliary[tag['TAG']] = tag
                                    else:
                                        log.error('ignoring invalid %s optional field %s', tag['TYPE'], tag['VALUE'])
                            segment['AUX'] = auxiliary

                        if self.read:
                            if self.read[0]['QNAME'] == segment['QNAME']:
                                self.read.append(segment)
                            else:
                                self.input_buffer.append(self.read)
                                self.read = [ segment ]
                                if len(self.input_buffer) >= self.buffer_capacity:
                                    break
                        else:
                            self.read = [ segment ]
                    else:
                        self.log.error('invalid sam syntax %s', line)

        if len(self.input_buffer) == 0 and self.read:
            self.input_buffer.append(self.read)
            self.read = None

        return len(self.input_buffer)

    def pick_barcode(self, decoder):
        barcode = None
        event = numpy.random.random()
        for record in decoder['distribution']:
            if event < record['density']:
                barcode = record['item']
                break
        if barcode:
            print('picked barcode {}'.format(barcode['index']))
        else:
            print('picked noise')
        return barcode

    def error_occured(self, position, nucleotide):
        if nucleotide != 'N':
            event = numpy.random.random() / self.model['error occurance factor']
            return event < self.model['error frequency'][position][nucleotide]
        else:
            return False

    def substitute_with(self, position, nucleotide):
        event = numpy.random.random()
        for option in self.model['substitution frequency'][position][nucleotide]:
            if event < option['density']:
                return option['to']
        return 'N'

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            for decoder in self.preset.values():
                barcode = self.pick_barcode(decoder)
                if barcode:
                    offset = 0
                    for token,barcode_segment in zip(decoder['transform']['parsed token'], barcode['barcode']):
                        read[token['segment']]['SEQ'][token['start']:token['end']] = barcode_segment
                        for p,i in zip(range(token['start'], token['end']), range(offset + token['length'])):
                            o = read[token['segment']]['SEQ'][p]
                            if self.error_occured(i, o):
                                s = self.substitute_with(i, o)
                                read[token['segment']]['SEQ'][p] = s
                                # print('substitute {} at position {} with {}'.format(o, p, s))
                        offset += token['length']

                    for segment in read:
                        segment['AUX'][decoder['true barcode tag']] = {
                            'TAG': decoder['true barcode tag'],
                            'TYPE': 'Z',
                            'VALUE': '-'.join([''.join(s) for s in barcode['barcode']])
                        }
                else:
                    # noise: insert a random piece of PhiX
                    pass

            self.output_buffer.append(read)

    def flush(self):
        while(len(self.output_buffer) > 0):
            read = self.output_buffer.pop(0)
            # print(read[0]['QNAME'])
            # print(to_json(read))
            for segment in read:
                line = '\t'.join([
                    segment['QNAME'],
                    str(segment['FLAG']),
                    segment['RNAME'],
                    str(segment['POS']),
                    str(segment['MAPQ']),
                    segment['CIGAR'],
                    segment['RNEXT'],
                    str(segment['PNEXT']),
                    str(segment['TLEN']),
                    ''.join(segment['SEQ']),
                    segment['QUAL']
                    # ''.join(segment['QUAL'])
                ])
                if 'AUX' in segment:
                    auxiliary = []
                    for tag in segment['AUX'].values():
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
                    line = '{}\t{}\n'.format(line, auxiliary)
                # self.output_handle.write(line)

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)

    pipeline = None
    try:
        pipeline = Simulate()
        pipeline.execute()

    except DownloadError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except ValueError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except BadConfigurationError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except CommandFailedError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except(KeyboardInterrupt, SystemExit) as e:
        if e.code != 0:
            logging.getLogger('main').critical(e)
            sys.exit(1)

    finally:
        if pipeline:
            pipeline.close()

    sys.exit(0)

if __name__ == '__main__':
    main()
