#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
from math import sqrt

# r['CC'],    # 0
# r['DX'],    # 1
# r['WL'],    # 2
# r['WO'],    # 3
# r['OG'],    # 4
# r['AS'],    # 5
# r['MQ'],    # 6
# r['QS'],    # 7
# r['OBS'],   # 8
# r['OBQ'],   # 9
# r['BEE'],   # 10
# r['PBS'],   # 11
# r['PBD'],   # 12
# r['DBS'],   # 13
# r['DBD'],   # 14
# r['P'],     # 15
# r['CP'],    # 16
# r['REE'],   # 17
# r['RS'],    # 18
# r['RQ'],    # 19
# r['ID'],    # 20
# r['BS'],    # 22
# r['BD'],    # 23

organism = [
    "BC",
    "CO",     # Y
    "NC",     # Y
    "PH",
    "PX",
    "SC",     # Y
    "SP",     # Y
    # "UK",
    "VR",
    # "Y"
]

MIN_KNOWN = 1000
MIN_YEAST = 100
KNOW_PORTION_MIN = 0.1
YEAST_PORTION_MIN = 0.1

EXCHANGE_MIN_KNOWN = 10
EXCHANGE_MIN_YEAST = 1
EXCHANGE_KNOW_PORTION_MIN = 0
EXCHANGE_YEAST_PORTION_MIN = 0


def to_json(node):
    print(json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4))

def pvalue():
    collection = json.load(sys.stdin)

    # annotate
    collection['order'] = sorted(list(collection['barcode'].keys()))
    for b in collection['order']:
        collection['barcode'][b]['key'] = b
        collection['barcode'][b]['error'] = {}
    collection['size'] = len(collection['order'])

    # adjust portion to not include UK
    for barcode in collection['barcode'].values():
        barcode['known count'] = 0
        for exchange in barcode['exchange'].values():
            exchange['known count'] = 0
            exchange['known portion'] = 0
            for k,v in exchange['organism'].items():
                if k != 'UK' and k != 'PX':
                    exchange['known count'] += v['count']

            for k,v in exchange['organism'].items():
                if k != 'UK':
                    v['portion'] = v['count'] / exchange['known count']

            if exchange['known count']:
                exchange['known portion'] = exchange['known count'] / exchange['known count']

            barcode['known count'] += exchange['known count']
            barcode['known portion'] = barcode['known count'] / barcode['count']

    # sum up 4 yeast samples into Y
    for barcode in collection['barcode'].values():
        barcode['yeast count'] = 0
        barcode['yeast portion'] = 0
        for exchange in barcode['exchange'].values():
            exchange['yeast count'] = 0
            exchange['yeast portion'] = 0
            for o in ['CO', 'NC', 'SC', 'SP']:
                if o in exchange['organism']:
                    exchange['yeast count'] += exchange['organism'][o]['count']
                    barcode['yeast count'] += exchange['organism'][o]['count']
            if exchange['known count']:
                exchange['yeast portion'] = exchange['yeast count'] / exchange['known count']
        if barcode['known count']:
            barcode['yeast portion'] = barcode['yeast count'] / barcode['known count']

    # compute self and other error for each exchange class
    for direction in ['1', '2']:

        r = 0
        while r < collection['size'] :
            right = collection['barcode'][collection['order'][r]]
            right['error'][direction] = { 'self': None, 'other':[] }

            l = 0
            while l < collection['size'] :
                left = collection['barcode'][collection['order'][l]]
                exchange_yeast_count = 0
                exchange_yeast_portion = 0
                exchange_known_count = 0
                exchange_known_portion = 0

                unchanged_yeast_count = 0
                unchanged_yeast_portion = 0
                unchanged_known_count = 0
                unchanged_known_portion = 0

                rss = 0

                exchange = None
                if direction in right['exchange']:
                    exchange = right['exchange'][direction]
                    exchange_yeast_count = exchange['yeast count']
                    exchange_yeast_portion = exchange['yeast portion']
                    exchange_known_count = exchange['known count']
                    exchange_known_portion = exchange['known portion']

                unchanged = None
                if '0' in left['exchange']:
                    unchanged = left['exchange']['0']
                    unchanged_yeast_count = unchanged['yeast count']
                    unchanged_yeast_portion = unchanged['yeast portion']
                    unchanged_known_count = unchanged['known count']
                    unchanged_known_portion = unchanged['known portion']

                if exchange and unchanged:
                    for o in organism:
                        portion_unchanged = 0
                        if o in unchanged['organism']:
                            portion_unchanged = unchanged['organism'][o]['portion']

                        portion_exchanged = 0
                        if o in exchange['organism']:
                            portion_exchanged = exchange['organism'][o]['portion']

                        rss += ((portion_unchanged - portion_exchanged)**2)

                elif unchanged:
                    for o in organism:
                        portion_unchanged = 0
                        if o in unchanged['organism']:
                            portion_unchanged = unchanged['organism'][o]['portion']
                        rss += (portion_unchanged**2)

                elif exchange:
                    for o in organism:
                        portion_exchanged = 0
                        if o in exchange['organism']:
                            portion_exchanged = exchange['organism'][o]['portion']
                        rss += (portion_exchanged**2)

                if l == r:
                    right['error'][direction]['self'] = sqrt(rss)

                elif (
                    exchange_yeast_count >= EXCHANGE_MIN_YEAST and
                    exchange_yeast_portion >= EXCHANGE_YEAST_PORTION_MIN and
                    exchange_known_count >= EXCHANGE_MIN_KNOWN and 
                    exchange_known_portion >= EXCHANGE_KNOW_PORTION_MIN and
                    unchanged_yeast_count >= EXCHANGE_MIN_YEAST and
                    unchanged_yeast_portion >= EXCHANGE_YEAST_PORTION_MIN and
                    unchanged_known_count >= EXCHANGE_MIN_KNOWN and 
                    unchanged_known_portion >= EXCHANGE_KNOW_PORTION_MIN and
                    left['known count'] >= MIN_KNOWN and
                    right['known count'] >= MIN_KNOWN and
                    left['known portion'] >= KNOW_PORTION_MIN and
                    right['known portion'] >= KNOW_PORTION_MIN and
                    left['yeast count'] >= MIN_YEAST and
                    right['yeast count'] >= MIN_YEAST and
                    left['yeast portion'] >= YEAST_PORTION_MIN and
                    right['yeast portion'] >= YEAST_PORTION_MIN
                ):
                    right['error'][direction]['other'].append([left['key'], sqrt(rss)])

                l += 1

            r += 1

        # compute mean and sd
        # for barcode in collection['barcode'].values():
        #     if 'rss other' in barcode:
        #         count = len(barcode['error'][direction]['other'])
        #         mean = sum(barcode['error'][direction]['other']) / count
        #         sd = sqrt(sum([(mean - x)**2 for x in barcode['error'][direction]['other']]) / count)
        #         barcode['error'][direction]['mean'] = mean
        #         barcode['error'][direction]['sd'] = sd

    return collection

def collect():
    collection = { 'barcode': {} }
    for line in sys.stdin:
        if line:
            row = line.strip().split()
            DX = row[1]
            OG = row[4]
            PBS = row[11]
            DBS = row[13]
            DBD = row[14]
            BS = row[21]
            BD = row[22]

            if '=' not in PBS or '=' not in DBS:
                if False:
                    # disabled because analyze_sam.py already computed DX correctly
                    EQ = 1 if DBS != PBS else 0
                    DX = 0
                    if EQ:
                        # REMOVED: probabilistic moved the read from a library to undetermined
                        if '=' in PBS:
                            DX = 1
                            BS = DBS
                            BD = DBD

                        # RECOVERED: probabilistic moved the read from a undetermined to a library
                        elif '=' in DBS:
                            DX = 2

                        # RECLASSIFIED: probabilistic moved the read from a one library to another
                        else:
                            DX = 3

                node = collection

                if BS not in node['barcode']:
                    node['barcode'][BS] =  { 'exchange': {}, 'count': 0 }
                node['barcode'][BS]['count'] += 1

                if DX not in node['barcode'][BS]['exchange']:
                    node['barcode'][BS]['exchange'][DX] = { 'organism': {}, 'count': 0 }
                node['barcode'][BS]['exchange'][DX]['count'] += 1

                if OG not in node['barcode'][BS]['exchange'][DX]['organism']:
                    node['barcode'][BS]['exchange'][DX]['organism'][OG] = { 'count': 0 }
                node['barcode'][BS]['exchange'][DX]['organism'][OG]['count'] += 1

    for barcode in collection['barcode'].values():
        for exchange in barcode['exchange'].values():
            for organism in exchange['organism'].values():
                organism['portion'] = organism['count'] / exchange['count']
        
        if 2 in barcode['exchange'] and 0 in barcode['exchange']:
            barcode['rss'] = 0
            for k,v in barcode['exchange'][0]['organism'].items():
                portion_unchanged = v['portion']
                portion_recovered = 0
                if k in barcode['exchange'][2]['organism']:
                    portion_recovered = barcode['exchange'][2]['organism'][k]['portion']
                barcode['rss'] += (portion_unchanged - portion_recovered)**2

    return collection

def csv(collection):
    buffer = []
    head = [
        'index',
        'DX',
        'base',
        'other',
        'self',
        'R',
    ]

    for direction in ['1', '2']:
        for index, key in enumerate(collection['order']):
            barcode = collection['barcode'][key]
            if direction in barcode['error']:
                row = [
                    index,
                    direction,
                    barcode['key'],
                    barcode['key'],
                    '1',
                    barcode['error'][direction]['self'],
                ]
                buffer.append(row)

                for value in barcode['error'][direction]['other']:
                    row = [
                        index,
                        direction,
                        barcode['key'],
                        value[0],
                        '0',
                        value[1],
                    ]
                    buffer.append(row)

    print('\t'.join([str(x) for x in head]))
    buffer.sort(key=lambda x: x[1], reverse=False)
    buffer.sort(key=lambda x: x[4], reverse=False)
    print('\n'.join(['\t'.join([str(x) for x in row]) for row in buffer]))

def phix(collection):
    phix_distribution = {
        'barcode': {},
        'phix total': 0,
        'phix total portion': 0,
        'phix mean': 0,
        'phix sd': 0,
    }
    for index, key in enumerate(collection['order']):
        barcode = collection['barcode'][key]
        phix_count = 0
        phix_portion = 0
        if 'PX' in barcode['exchange']['0']['organism']:
            phix_count = barcode['exchange']['0']['organism']['PX']['count']
            phix_portion = barcode['exchange']['0']['organism']['PX']['portion']
            phix_distribution['phix total'] += phix_count
            phix_distribution['phix total portion'] += phix_portion

        phix_distribution['barcode'][key] = {
            'barcode': key,
            'count': barcode['count'],
            'phix count': phix_count,
            'phix portion': phix_portion,
        }

    count = len(phix_distribution['barcode'])
    mean = phix_distribution['phix total portion'] / len(phix_distribution['barcode'])
    sd = sqrt(sum([(mean - x['phix portion'])**2 for x in phix_distribution['barcode'].values()]) / count)
    phix_distribution['phix mean'] = mean
    phix_distribution['phix sd'] = sd
    to_json(phix_distribution)
    # print('\t'.join([str(x) for x in [barcode['key'], phix_count, phix_portion]]))


# collection = collect()
collection = pvalue()
csv(collection)
# to_json(collection)
# phix(collection)