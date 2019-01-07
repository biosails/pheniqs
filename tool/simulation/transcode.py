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

import io
import os
import sys
import json
import numpy
from copy import deepcopy
from threading import Thread
from queue import Queue, Empty
from subprocess import Popen, PIPE

from core.error import *
from core import Job
from core import merge
from core import to_json
from core import remove_compiled
from core import prepare_path


class BamReader(object):
    def __init__(self, path):
        self.eof = False
        self.queue = Queue(128)
        self.input_process = Popen(args=[ 'samtools', 'view', path ], stdout=PIPE, stderr=PIPE)
        self.thread = Thread(target = self.replenish)
        self.thread.start()

    def replenish(self):
        for line in iter(self.input_process.stdout.readline, b''):
            if line:
                self.queue.put(line.strip().decode('utf-8'))
            else:
                self.eof = True

    def read_line(self):
        line = None
        while line is None and not self.eof:
            try:
                line = self.queue.get(timeout = 0.1)
            except Empty:
                break
        return line

class BamWriter(object):
    def __init__(self, path):
        self.eof = False
        self.queue = Queue(128)
        self.output_process = Popen(args=[ 'samtools', 'view', '-b', '-o', path], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        self.thread = Thread(target = self.flush)
        self.thread.start()

    def flush(self):
        while True:
            try:
                line = self.queue.get(timeout = 0.1)
                self.output_process.stdin.write(line.encode('utf-8'))
                self.output_process.stdin.write(b'\n')
            except Empty:
                if self.eof: break

    def write_line(self, line):
        self.queue.put(line)

    def close(self):
        self.eof = True

class Transcode(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        default = {
            'phred scale': [],
            'decoder by index': [],
        }
        self.ontology = merge(default, self.ontology)
        self.input_buffer = None
        self.input_feed = None
        self.output_buffer = None
        self.output_feed = None
        self.read_buffer = None

    @property
    def buffer_capacity(self):
        return self.instruction['capacity']

    @property
    def decoder_by_index(self):
        return self.ontology['decoder by index']

    @property
    def phred_scale(self):
        return self.ontology['phred scale']

    @property
    def instrument_id(self):
        return self.model['instrument id']

    @property
    def flowcell_id(self):
        return self.model['flowcell id']

    @property
    def run_count(self):
        return self.model['run count']

    def load_phred_scale(self):
        phred_probability_base = pow(10.0, -0.1)
        for q in range(64):
            self.ontology['phred scale'].append({
                'code': chr(q + 33),
                'value': q,
                'error': pow(phred_probability_base, q)
            })

    def parse_qname(self, read):
        # @PHENIQS:1:HABDFADXX:0112241932:23:31:4
        #<instrument id>:<run count>:<flowcell id>:<read index>[:<decoder hint>]*
        read['hint'] = []
        parsed = read['segment'][0]['fixed'][0].split(':')

        if len(parsed) - 4 == len(self.decoder_by_index):
            for field in parsed[4:]:
                read['hint'].append(int(field))
        else:
            raise SequenceError('incorrect number of decoder hints')

    def load(self):
        self.load_phred_scale()

    def open(self):
        if os.path.exists(self.instruction['input']):
            self.input_feed = BamReader(self.instruction['input'])
            self.input_buffer = []
        else:
            raise BadConfigurationError('input {} not found'.format(self.instruction['input']))

        if 'output' in self.instruction:
            prepare_path(os.path.abspath(os.path.expanduser(os.path.expandvars(self.instruction['output']))), self.log)
            self.output_feed = BamWriter(self.instruction['output'])
            self.output_buffer = []

    def close(self):
        Job.close(self)
        if self.output_feed:
            self.output_feed.close()

    def execute(self):
        self.load()
        self.open()
        while(self.replenish()):
            self.manipulate()
            self.flush()
        self.close()
        self.finalize()
        if self.output_feed:
            self.output_feed.thread.join()

    def replenish(self):
        while True:
            line = self.input_feed.read_line()
            if line is not None:
                if line[0] =='@':
                    # redirect header line
                    self.output_feed.write_line(line)
                else:
                    segment = { 'fixed': line.split('\t'), 'auxiliary': {} }
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
            else:
                break

        if self.input_feed.eof and len(self.input_buffer) == 0 and self.read_buffer:
            self.input_buffer.append(self.read_buffer)
            self.read_buffer = None

        return len(self.input_buffer)

        # for line in self.input_feed:
        #     if line:
        #         if line[0] =='@':
        #             # redirect header line
        #             self.output_feed.write(line)
        #         else:
        #             segment = { 'fixed': line.strip().split('\t'), 'auxiliary': {} }
        #             if len(segment['fixed']) > 10:
        #                 # parse AUX
        #                 if len(segment['fixed']) > 11:
        #                     for item in segment['fixed'][11:]:
        #                         if len(item) > 5:
        #                             tag = { 'TAG': item[0:2], 'TYPE': item[3:4], 'VALUE': item[5:] }
        #                             if tag['TYPE'] is 'i':
        #                                 tag['VALUE'] = int(tag['VALUE'])
        #                             elif tag['TYPE'] is 'f':
        #                                 tag['VALUE'] = float(tag['VALUE'])
        #                             segment['auxiliary'][tag['TAG']] = tag
        #                         else:
        #                             self.log.warning('ignoring invalid auxiliary tag %s', field)
        #
        #                 if self.read_buffer:
        #                     if self.read_buffer['segment'][0]['fixed'][0] == segment['fixed'][0]:
        #                         # if QNAME is the same as the first segment in the read buffer
        #                         # this is another segment of the same read
        #                         self.read_buffer['segment'].append(segment)
        #                     else:
        #                         # otherwise this is a segment from a new read
        #                         # so write the completed read to the input buffer and intitalize the read buffer
        #                         self.input_buffer.append(self.read_buffer)
        #                         self.read_buffer = { 'segment': [ segment ] }
        #                 else:
        #                     # only when the first read segment is encountered
        #                     self.read_buffer = { 'segment': [ segment ] }
        #             else:
        #                 # a sam record must have at least 11 mandatory fields
        #                 self.log.error('invalid sam syntax %s', line)
        #                 raise SequenceError(line)
        #
        #             if len(self.input_buffer) >= self.buffer_capacity: break
        #
        # if len(self.input_buffer) == 0 and self.read_buffer:
        #     self.input_buffer.append(self.read_buffer)
        #     self.read_buffer = None
        #
        # return len(self.input_buffer)

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
                self.output_feed.write_line('{}\t{}'.format(fixed, auxiliary))

    def finalize(self):
        pass

    def parse_rule(self, rule):
        node = {}
        if 'token' in rule:
            node['token'] = []
            for pattern in rule['token']:
                parsed = pattern.split(':')
                token = { 'segment': parsed[0], 'start': parsed[1], 'end': parsed[2] }

                if len(token['segment']) > 0:
                    token['segment'] = int(token['segment'])
                else:
                    raise BadConfigurationError('segment index is missing in {}'.format(pattern))

                if len(token['start']) > 0:
                    token['start'] = int(token['start'])
                else:
                    token['start'] = 0

                if len(token['end']) > 0:
                    token['end'] = int(token['end'])
                else:
                    raise BadConfigurationError('end position is missing in {}'.format(pattern))

                token['length'] = token['end'] - token['start']
                node['token'].append(token)
        return node
