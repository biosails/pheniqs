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

class SamTranscode(Pipeline):
    def __init__(self, name):
        Pipeline.__init__(self, name)
        self.read_buffer = None
        self.input_buffer = None
        self.output_buffer = None
        self.input_handle = None
        self.output_handle = None
        self.ontology['token model'] = {}
        self.ontology['decoder by index'] = []

    @property
    def preset(self):
        return self.ontology['preset']

    @property
    def buffer_capacity(self):
        return self.instruction['capacity']

    @property
    def instrument_id(self):
        return self.preset['instrument id']

    @property
    def flowcell_id(self):
        return self.preset['flowcell id']

    @property
    def run_count(self):
        return self.preset['run count']

    @property
    def read_count(self):
        return self.preset['count']

    @property
    def token_model(self):
        return self.ontology['token model']

    def decoder_by_index(self):
        return self.ontology['decoder by index']

    def load(self):
        pass

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

    def execute(self):
        self.load()
        self.open()
        while(self.replenish()):
            self.manipulate()
            self.flush()
        self.close()

        self.finalize()
        self.report()

    def replenish(self):
        for line in self.input_handle:
            if line:
                if line[0] =='@':
                    # redirect header line
                    self.output_handle.write(line)
                else:
                    segment = { 'fixed': line.strip().split('\t'), 'auxiliary': {} }
                    if len(segment['fixed']) > 10:
                        # parse AUX
                        if len(segment['fixed']) > 11:
                            for item in segment['fixed'][11:]:
                                if len(item) > 5:
                                    tag = { 'TAG': item[0:2], 'TYPE': item[3:4], 'VALUE': item[5:] }
                                    if tag['TYPE'] is 'i':
                                        tag['VALUE'] = int(tag['VALUE'])
                                    elif tag['TYPE'] is 'f':
                                        tag['VALUE'] = float(tag['VALUE'])
                                    segment['auxiliary'][tag['TAG']] = tag
                                else:
                                    self.log.warning('ignoring invalid auxiliary tag %s', field)

                        # segment['fixed'][9] = numpy.array(list(segment['fixed'][9]))
                        if self.read_buffer:
                            if self.read_buffer['segment'][0]['fixed'][0] == segment['fixed'][0]:
                                # if QNAME is the same as the first segment in the read buffer
                                # this is another segment of the same read
                                self.read_buffer['segment'].append(segment)
                            else:
                                # otherwise this is a segment from a new read
                                # so write the completed read to the input buffer and intitalize the read buffer
                                self.input_buffer.append(self.read_buffer)
                                self.read_buffer = { 'segment': [ segment ] }
                        else:
                            # only when the first read segment is encountered
                            self.read_buffer = { 'segment': [ segment ] }
                    else:
                        # a sam record must have at least 11 mandatory fields
                        self.log.error('invalid sam syntax %s', line)
                        raise SequenceError(line)

                    if len(self.input_buffer) >= self.buffer_capacity: break

        if len(self.input_buffer) == 0 and self.read_buffer:
            self.input_buffer.append(self.read_buffer)
            self.read_buffer = None

        return len(self.input_buffer)

    def manipulate(self):
        pass

    def flush(self):
        while(len(self.output_buffer) > 0):
            read = self.output_buffer.pop(0)
            for segment in read['segment']:
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
                    segment['fixed'][9],
                    segment['fixed'][10]
                ])
                if segment['auxiliary']:
                    auxiliary = []
                    for tag in segment['auxiliary'].values():
                        if tag['TYPE'] is 'A':
                            pass
                        elif tag['TYPE'] is 'i':
                            tag['VALUE'] = str(tag['VALUE'])
                        elif tag['TYPE'] is 'f':
                            tag['VALUE'] = str(tag['VALUE'])
                        elif tag['TYPE'] is 'Z':
                            pass
                        elif tag['TYPE'] is 'H':
                            pass
                        elif tag['TYPE'] is 'B':
                            pass
                        auxiliary.append('{}:{}:{}'.format(tag['TAG'], tag['TYPE'], tag['VALUE']))
                    auxiliary = '\t'.join(auxiliary)
                self.output_handle.write('{}\t{}\n'.format(fixed, auxiliary))

    def finalize(self):
        pass

    def report(self):
        print_json(self.preset)

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
