#!/usr/bin/env python3

# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2017  Lior Galanti
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

import io
import json
import sys
import re
from datetime import timedelta

def to_json(node):
    print(json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4))

sample_expression = re.compile(r'(?P<user>[0-9\.]+)s user (?P<system>[0-9\.]+)s system (?P<cpu>[0-9]+)% cpu (?P<total>[0-9\.:]+) total')
expression = re.compile(r'^(?P<flowcell>HCJFNBCXX|HG7CVAFXX|HK5NHBGXX)\sbenchmark\s(?P<utility>pheniqs|fastq-multx|picard)\s(?P<decoder>mdd|pamld)\s(?P<layout>split|interleaved|combined)\s(?P<format>fastq|cram)(?:\s(?P<quality>quality))?$')

batch = {}
node = {}
key = None

for line in sys.stdin:
    line = line.strip()
    if line and line[0] != '#':
        match = expression.search(line)
        if match:
            key = line;
            if key not in batch:
                record = match.groupdict()
                record['key'] = line
                record['sample'] = []
                if 'quality' in record and record['quality']:
                    record['quality'] = True
                else:
                    record['quality'] = False
                batch[key] = record
        elif key is not None:
            match = sample_expression.search(line)
            if match:
                batch[key]['sample'].append(match.groupdict())
        else:
            print('encountered unknown line {}'.format(line))
            key = None
            break

for key, value in batch.items():
    value['sample count'] = 0
    value['cpu'] = 0
    value['system'] = 0
    value['total'] = 0
    value['user'] = 0
    for sample in value['sample']:
        value['sample count'] += 1
        value['cpu'] += float(sample['cpu'])
        value['system'] += float(sample['system'])
        # value['total'] += sample['total'] # 6:04.95
        value['user'] += float(sample['user'])

    value['cpu'] /= value['sample count']
    value['system'] /= value['sample count']
    # value['total'] /= value['sample count']
    value['user'] /= value['sample count']
    value['total'] = round((value['user'] + value['system']) / value['cpu'] * 100 )
    value['total time'] = str(timedelta(seconds=value['total']))

for record in batch.values():
    pivot = node
    if record['flowcell'] not in pivot:
        pivot[record['flowcell']] = {}

    pivot = pivot[record['flowcell']]
    if record['decoder'] not in pivot:
        pivot[record['decoder']] = {}

    pivot = pivot[record['decoder']]
    if record['format'] not in pivot:
        pivot[record['format']] = {}

    pivot = pivot[record['format']]
    if record['layout'] not in pivot:
        pivot[record['layout']] = {}

    pivot = pivot[record['layout']]
    if record['utility'] not in pivot:
        pivot[record['utility']] = { 
            'with quality': None,
            'without quality': None
        }
    pivot = pivot[record['utility']]
    if record['quality']:
        pivot['with quality'] = record
    else:
        pivot['without quality'] = record


def print_R_table(node):
    print('\t'.join([
        'description',
        'flowcell',
        'utility',
        'decoder',
        'layout',
        'format',
        'seconds',
        'cpu',
        'system',
        'user',
        'quality'
    ]))
    for flowcell, flowcell_node in node.items():
        for decoder, decoder_node in flowcell_node.items():
            for format, format_node in decoder_node.items():
                for layout, layout_node in format_node.items():
                    for utility, pivot in layout_node.items():
                        nq = pivot['without quality']
                        q = pivot['with quality']
                        pivot['total'] = nq['total']
                        pivot['cpu'] = nq['cpu']
                        pivot['system'] = nq['system']
                        pivot['user'] = nq['user']
                        pivot['utility'] = nq['utility']
                        pivot['flowcell'] = nq['flowcell']
                        pivot['utility'] = nq['utility']
                        pivot['decoder'] = nq['decoder']
                        pivot['layout'] = nq['layout']
                        pivot['format'] = nq['format']
                        pivot['delta total'] = 0
                        pivot['delta cpu'] = 0
                        pivot['delta system'] = 0
                        pivot['delta user'] = 0
                        pivot['description'] = ' '.join ([
                            pivot['utility'],
                            pivot['decoder'],
                            pivot['layout'],
                            pivot['format'],
                        ])

                        if q:
                            pivot['delta total'] = max(0, q['total'] - nq['total'])
                            pivot['delta cpu'] = max(0, q['cpu'] - nq['cpu'])
                            pivot['delta system'] = max(0, q['system'] - nq['system'])
                            pivot['delta user'] = max(0, q['user'] - nq['user'])

                        record = [
                            pivot['description'],
                            pivot['flowcell'],
                            pivot['utility'],
                            pivot['decoder'],
                            pivot['layout'],
                            pivot['format'],
                            pivot['total'],
                            pivot['cpu'],
                            pivot['system'],
                            pivot['user'],
                            1
                        ]

                        record = [ str(v) for v in record ]
                        print('\t'.join(record))

                        record = [
                            pivot['description'],
                            pivot['flowcell'],
                            pivot['utility'],
                            pivot['decoder'],
                            pivot['layout'],
                            pivot['format'],
                            pivot['delta total'],
                            pivot['delta cpu'],
                            pivot['delta system'],
                            pivot['delta user'],
                            0
                        ]

                        record = [ str(v) for v in record ]
                        print('\t'.join(record))

def print_markdown_table(node):
    print('\t'.join([
        'description',
        'flowcell',
        'cpu',
        'seconds',
        'q delta',
    ]))
    for flowcell, flowcell_node in node.items():
        for decoder, decoder_node in flowcell_node.items():
            for format, format_node in decoder_node.items():
                for layout, layout_node in format_node.items():
                    for utility, pivot in layout_node.items():
                        nq = pivot['without quality']
                        q = pivot['with quality']
                        pivot['total'] = nq['total']
                        pivot['cpu'] = nq['cpu']
                        pivot['system'] = nq['system']
                        pivot['user'] = nq['user']
                        pivot['utility'] = nq['utility']
                        pivot['flowcell'] = nq['flowcell']
                        pivot['utility'] = nq['utility']
                        pivot['decoder'] = nq['decoder']
                        pivot['layout'] = nq['layout']
                        pivot['format'] = nq['format']
                        pivot['delta total'] = 0
                        pivot['delta cpu'] = 0
                        pivot['delta system'] = 0
                        pivot['delta user'] = 0
                        pivot['description'] = ' '.join ([
                            pivot['utility'],
                            pivot['decoder'],
                            pivot['layout'],
                            pivot['format'],
                        ])

                        if q:
                            pivot['delta total'] = max(0, q['total'] - nq['total'])
                            pivot['delta cpu'] = max(0, q['cpu'] - nq['cpu'])
                            pivot['delta system'] = max(0, q['system'] - nq['system'])
                            pivot['delta user'] = max(0, q['user'] - nq['user'])

                        if(pivot['cpu'] > 0):
                            record = [
                                pivot['description'],
                                pivot['flowcell'],
                                int(pivot['cpu']),
                                timedelta(seconds=pivot['total']),
                                int(pivot['delta total']),
                            ]

                            record = [ str(v) for v in record ]
                            print('|'.join(record))

                            # record = [
                            #     pivot['description'],
                            #     pivot['flowcell'],
                            #     int(pivot['delta total']),
                            #     int(pivot['delta cpu']),
                            #     int(pivot['delta system']),
                            #     int(pivot['delta user']),
                            #     pivot['delta system']/pivot['delta user'],
                            # ]

                            # record = [ str(v) for v in record ]
                            # print('|'.join(record))

def print_latex_table(node):
    print('\t'.join([
        'Description',
        'Layout',
        'Seconds',
        'Delta',
    ]))
    for flowcell, flowcell_node in node.items():
        for decoder, decoder_node in flowcell_node.items():
            for format, format_node in decoder_node.items():
                for layout, layout_node in format_node.items():
                    for utility, pivot in layout_node.items():
                        nq = pivot['without quality']
                        q = pivot['with quality']
                        pivot['total'] = nq['total']
                        pivot['cpu'] = nq['cpu']
                        pivot['system'] = nq['system']
                        pivot['user'] = nq['user']
                        pivot['utility'] = nq['utility']
                        pivot['flowcell'] = nq['flowcell']
                        pivot['utility'] = nq['utility']
                        pivot['decoder'] = nq['decoder']
                        pivot['layout'] = nq['layout']
                        pivot['format'] = nq['format']
                        pivot['delta total'] = 0
                        pivot['delta cpu'] = 0
                        pivot['delta system'] = 0
                        pivot['delta user'] = 0
                        pivot['description'] = ' '.join ([
                            pivot['utility'],
                            pivot['decoder'],
                            pivot['layout'],
                            pivot['format'],
                        ])

                        if(pivot['flowcell'] == 'HCJFNBCXX'):
                            pivot['flowcell layout'] = '2x250PE'
                        elif(pivot['flowcell'] == 'HG7CVAFXX'):
                            pivot['flowcell layout'] = '2x75PE'
                        elif(pivot['flowcell'] == 'HK5NHBGXX'):
                            pivot['flowcell layout'] = '2x36PE'

                        if q:
                            pivot['delta total'] = max(0, q['total'] - nq['total'])
                            pivot['delta cpu'] = max(0, q['cpu'] - nq['cpu'])
                            pivot['delta system'] = max(0, q['system'] - nq['system'])
                            pivot['delta user'] = max(0, q['user'] - nq['user'])

                        if(pivot['cpu'] > 0):
                            record = [
                                pivot['description'],
                                pivot['flowcell layout'],
                                timedelta(seconds=pivot['total']),
                                timedelta(seconds=pivot['delta total']),
                            ]

                            record = [ str(v) for v in record ]
                            print('{}\t{}'.format('\t & '.join(record), r'\\'))

print_R_table(node)
# print_latex_table(node)
# print_markdown_table(node)
# to_json(node)

    # print(','.join(['utility', 'decoder', 'layout', 'format', 'quality', 'cpu', 'system', 'user', 'time']))
    # for key, value in batch.items():
    #     print (
    #         ','.join (
    #             [
    #                 value['utility'],
    #                 value['decoder'],
    #                 value['layout'],
    #                 value['format'],
    #                 str(value['quality']),
    #                 '{:10.0f}'.format(value['cpu']),
    #                 '{:10.0f}'.format(value['system']),
    #                 '{:10.0f}'.format(value['user']),
    #                 value['total'].lstrip('0:')
    #             ]
    #         )
    #     )
