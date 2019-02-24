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

import logging
from core.error import *
from core import Job
from core import to_json

def fscore(precision, recall):
    if(precision + recall) > 0:
        return 2.0 * (precision * recall) / (precision + recall)
    else:
        return 0

class Summarize(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.log = logging.getLogger('Summarize')

    @property
    def experiment(self):
        return self.ontology['experiment']

    @property
    def collection(self):
        return self.experiment['collection']

    def collect_accuracy_benchmark(self):
        # Noise accumulator:
        #   TP :
        #       noise: a noise read was correctly classified as noise
        #   FN :
        #       noise: a noise read was incorrectly classified to a real barcode
        #   FP :
        #       noise: a real read was classified as noise
        #
        # Real barcode accumulator:
        #   TP :
        #       real: real read is correctly classified
        #   FN :
        #       real: a read from this barcode was incorrectly classified to another real barcode
        #       noise: a read from this barcode was incorrectly classified as noise
        #   FP :
        #       real: a real read from another barcode was incorrectly classified to this barcode
        #       noise: a noise read was incorrectly classified to this barcode
        def collect_barcode(n, barcode):
            d = n['decoder']
            b = n['barcode'][barcode['index']]
            for rank in [ 'noise', 'real' , 'both' ]:
                b[rank] = {}
                for qc in ['fail', 'pass', 'both']:
                    b[rank][qc] = {}
                    for item in [ 'count', 'TP', 'FP', 'FN', 'FDR', 'MR', 'precision', 'recall', 'fscore' ]:
                        b[rank][qc][item] = 0

            for rank in [ 'noise', 'real' ]:
                for qc in ['fail', 'pass']:
                    for item in [ 'TP', 'FP', 'FN' ]:
                        b[rank][qc][item] = barcode['accumulate'][rank][qc][item]
                        b[rank]['both'][item] += b[rank][qc][item]
                        b['both'][qc][item] += barcode['accumulate'][rank][qc][item]
                        b['both']['both'][item] += b[rank][qc][item]

                        if b['index'] > 0:
                            d[rank][qc][item] += b[rank][qc][item]
                            d[rank]['both'][item] += b[rank][qc][item]
                            d['both'][qc][item] += b[rank][qc][item]
                            d['both']['both'][item] += b[rank][qc][item]

                    b[rank][qc]['count'] = b[rank][qc]['TP'] + b[rank][qc]['FN']
                    b[rank]['both']['count'] += b[rank][qc]['count']
                    b['both'][qc]['count'] += b[rank][qc]['count']
                    b['both']['both']['count'] += b[rank][qc]['count']

                    if b['index'] > 0:
                        d[rank][qc]['count'] += b[rank][qc]['count']
                        d[rank]['both']['count'] += b[rank][qc]['count']
                        d['both'][qc]['count'] += b[rank][qc]['count']
                        d['both']['both']['count'] += b[rank][qc]['count']

            for rank in [ 'noise', 'real', 'both' ]:
                for qc in ['fail', 'pass', 'both']:
                    if b[rank][qc]['TP'] > 0 or b[rank][qc]['FP'] > 0:
                        b[rank][qc]['precision']    = b[rank][qc]['TP'] / (b[rank][qc]['FP'] + b[rank][qc]['TP'])
                        b[rank][qc]['FDR']          = b[rank][qc]['FP'] / (b[rank][qc]['FP'] + b[rank][qc]['TP'])
                    if b[rank][qc]['FN'] > 0 or b[rank][qc]['TP'] > 0:
                        b[rank][qc]['recall']       = b[rank][qc]['TP'] / (b[rank][qc]['TP'] + b[rank][qc]['FN'])
                        b[rank][qc]['MR']           = b[rank][qc]['FN'] / (b[rank][qc]['TP'] + b[rank][qc]['FN'])
                    b[rank][qc]['fscore'] = fscore(b[rank][qc]['precision'], b[rank][qc]['recall'])

        def collect_decoder(n, decoder):
            n['decoder'] = {}
            n['barcode'] = []
            for index in range(decoder['barcode cardinality'] + 1):
                b = { 'index': index }
                n['barcode'].append(b)

            d = n['decoder']
            for rank in [ 'noise', 'real', 'both' ]:
                d[rank] = {}
                for qc in ['fail', 'pass', 'both']:
                    d[rank][qc] = {}
                    for item in [ 'count', 'TP', 'FP', 'FN', 'FDR', 'MR', 'precision', 'recall', 'fscore' ]:
                        d[rank][qc][item] = 0

            collect_barcode(n, decoder['unclassified'])
            for barcode in decoder['codec'].values():
                collect_barcode(n, barcode)

            for rank in [ 'noise', 'real', 'both' ]:
                for qc in ['fail', 'pass', 'both']:
                    if d[rank][qc]['TP'] > 0 or d[rank][qc]['FP'] > 0:
                        d[rank][qc]['precision']    = d[rank][qc]['TP'] / (d[rank][qc]['FP'] + d[rank][qc]['TP'])
                        d[rank][qc]['FDR']          = d[rank][qc]['FP'] / (d[rank][qc]['FP'] + d[rank][qc]['TP'])
                    if d[rank][qc]['FN'] > 0 or d[rank][qc]['TP'] > 0:
                        d[rank][qc]['recall']       = d[rank][qc]['TP'] / (d[rank][qc]['TP'] + d[rank][qc]['FN'])
                        d[rank][qc]['MR']           = d[rank][qc]['FN'] / (d[rank][qc]['TP'] + d[rank][qc]['FN'])
                    d[rank][qc]['fscore'] = fscore(d[rank][qc]['precision'], d[rank][qc]['recall'])

        if 'collection' not in  self.experiment:
            self.log.info('computing collection')
            collection = []
            for substitution in self.experiment['substitution'].values():
                if 'model' in substitution:
                    record = {
                        'expected': substitution['model']['multiplex']['expected substitution rate'],
                        'simulated': substitution['model']['multiplex']['simulated substitution rate'],
                        'benchmark': []
                    }
                    if 'calibrated' in substitution['model'] and substitution['model']['calibrated']:
                        record['requested'] = substitution['model']['multiplex']['requested substitution rate'],


                    for tool, akey in {
                        'deml': 'deml demultiplex analysis',
                        'pamld': 'pamld demultiplex analysis',
                        'pamld_ap': 'pamld accurate prior demultiplex analysis',
                        'pamld_u': 'pamld uniform demultiplex analysis',
                        'mdd': 'mdd demultiplex analysis',
                    }.items():
                        if akey in substitution:
                            analysis = { 'tool': tool, }
                            if 'multiplex' in substitution[akey]:
                                if 'accumulate' in substitution[akey]['multiplex']['unclassified']:
                                    analysis['multiplex'] = {}
                                    collect_decoder(analysis['multiplex'], substitution[akey]['multiplex'])

                            if 'multiplex' in analysis:
                                record['benchmark'].append(analysis)

                    if record['benchmark']:
                        collection.append(record)

            collection.sort(key=lambda i: i['simulated'])
            self.experiment['collection'] = []

            current = 0
            for record in collection:
                test = abs(1.0 - current / record['simulated'])
                if test > 0.01:
                    current = record['simulated']
                    self.experiment['collection'].append(record)

    def summarize_decoder_accuracy_benchmark(self):
        header = [
            'rate',
            'tool',
            'rank',
            'qc',
            'TP',
            'FP',
            'FN', # for real ranl reads only FP and FN are the same
            # 'FDR',
            # 'MR',
            # 'precision',
            # 'recall',
        ]
        collection = []
        for record in self.collection:
            rate = record['simulated']
            if rate < 0.06:
                for benchmark in record['benchmark']:
                    tool = benchmark['tool']
                    for rank in [ 'noise', 'real', 'both' ]:
                        for qc in [ 'pass', 'fail', 'both' ]:
                        # for qc in [ 'pass', 'fail', 'both' ]:
                            qnode = benchmark['multiplex']['decoder'][rank][qc]
                            record = [
                                rate,
                                tool,
                                rank,
                                qc,
                                qnode['TP'],
                                qnode['FP'],
                                qnode['FN'],
                                # qnode['FDR'],
                                # qnode['MR'],
                                # qnode['precision'],
                                # qnode['recall'],
                            ]
                            collection.append(record)

        collection.sort(key=lambda i: i[3])
        collection.sort(key=lambda i: i[2])
        collection.sort(key=lambda i: i[1])
        collection.sort(key=lambda i: i[0])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def summarize_barcode_accuracy_benchmark(self):
        header = [
            'index',
            'rate',
            'tool',
            'rank',
            'qc',
            'TP',
            'FP',
            'FN',
            'FDR',
            'MR',
            # 'precision',
            # 'recall',
        ]
        collection = []
        for record in self.collection:
            rate = record['simulated']
            for benchmark in record['benchmark']:
                tool = benchmark['tool']
                for rank in [ 'noise', 'real' ]:
                # for rank in [ 'noise', 'real', 'both' ]:
                    for qc in [ 'pass', 'fail', 'both' ]:
                    # for qc in [ 'pass', 'fail', 'both' ]:
                        for barcode in benchmark['multiplex']['barcode']:
                            if barcode['index'] > 0:
                                qnode = barcode[rank][qc]
                                record = [
                                    barcode['index'],
                                    rate,
                                    tool,
                                    rank,
                                    qc,
                                    qnode['TP'],
                                    qnode['FP'],
                                    qnode['FN'],
                                    qnode['FDR'],
                                    qnode['MR'],
                                ]
                                collection.append(record)

        collection.sort(key=lambda i: i[4])
        collection.sort(key=lambda i: i[3])
        collection.sort(key=lambda i: i[2])
        collection.sort(key=lambda i: i[1])
        collection.sort(key=lambda i: i[0])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def summarize_noise_accuracy_benchmark(self):
        header = [
            'rate',
            'tool',
            'qc',
            'TP',
            'FP',
            'FN',
            'FDR',
            'MR',
            'precision',
            'recall',
        ]
        collection = []
        for record in self.collection:
            rate = record['simulated']
            for benchmark in record['benchmark']:
                tool = benchmark['tool']
                for qc in [ 'pass', 'fail', 'both' ]:
                # for qc in [ 'pass', 'fail', 'both' ]:
                    qnode = benchmark['multiplex']['barcode'][0]['noise'][qc]
                    record = [
                        rate,
                        tool,
                        qc,
                        qnode['TP'],
                        qnode['FP'],
                        qnode['FN'],
                        qnode['FDR'],
                        qnode['MR'],
                        qnode['precision'],
                        qnode['recall'],
                    ]
                    collection.append(record)

        collection.sort(key=lambda i: i[3])
        collection.sort(key=lambda i: i[2])
        collection.sort(key=lambda i: i[1])
        collection.sort(key=lambda i: i[0])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def summarize_decoder_quality_distribution_R(self):
        header = [
            'rate',
            'quality',
            'density',
        ]
        collection = []
        for substitution in self.experiment['substitution'].values():
            if 'model' in substitution:
                record = {
                    'expected': substitution['model']['multiplex']['expected substitution rate'],
                    'simulated': substitution['model']['multiplex']['simulated substitution rate'],
                    'benchmark': []
                }
                rate = record['simulated']
                if 'multiplex' in substitution['model']:
                    if 'quality distribution' in substitution['model']['multiplex']:
                        for quality, density in enumerate(substitution['model']['multiplex']['quality distribution']):
                            record = [
                                rate,
                                quality,
                                density
                            ]
                            collection.append(record)

        collection.sort(key=lambda i: i[1])
        collection.sort(key=lambda i: i[0])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def summarize_decoder_accuracy_benchmark_R(self):
        header = [
            'rate',
            'tool',
            'rank',
            'qc',
            'variable',
            'value',
        ]
        collection = []
        for record in self.collection:
            rate = record['simulated']
            for benchmark in record['benchmark']:
                tool = benchmark['tool']
                for rank in [ 'both' ]:
                # for rank in [ 'noise', 'real', 'both' ]:
                    for qc in [ 'pass', 'fail', 'both' ]:
                    # for qc in [ 'pass', 'fail', 'both' ]:
                        qnode = benchmark['multiplex']['decoder'][rank][qc]
                        for variable in [ 'FDR', 'MR', 'TP', 'FN', 'FP', 'fscore' ]:
                        # for variable in [ 'TP', 'FP', 'FN', 'FDR', 'MR', 'precision', 'recall' ]:
                            record = [
                                rate,
                                tool,
                                rank,
                                qc,
                                variable,
                                qnode[variable],
                            ]
                            collection.append(record)

        collection.sort(key=lambda i: i[3])
        collection.sort(key=lambda i: i[2])
        collection.sort(key=lambda i: i[1])
        collection.sort(key=lambda i: i[0])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def summarize_barcode_accuracy_benchmark_R(self):
        header = [
            'index',
            'rate',
            'tool',
            'rank',
            'qc',
            'variable',
            'value',
        ]
        collection = []
        for record in self.collection:
            rate = record['simulated']
            for benchmark in record['benchmark']:
                tool = benchmark['tool']
                for rank in [ 'real', 'noise' ]:
                # for rank in [ 'noise', 'real', 'both' ]:
                    for qc in [ 'pass', 'fail', 'both' ]:
                    # for qc in [ 'pass', 'fail', 'both' ]:
                        for barcode in benchmark['multiplex']['barcode']:
                            if barcode['index'] > 0:
                                qnode = barcode[rank][qc]
                                for variable in [ 'TP', 'FP', 'FN', 'FDR', 'MR' ]:
                                    record = [
                                        barcode['index'],
                                        rate,
                                        tool,
                                        rank,
                                        qc,
                                        variable,
                                        qnode[variable],
                                    ]
                                    collection.append(record)

        collection.sort(key=lambda i: i[4])
        collection.sort(key=lambda i: i[3])
        collection.sort(key=lambda i: i[2])
        collection.sort(key=lambda i: i[1])
        collection.sort(key=lambda i: i[0])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def summarize_noise_accuracy_benchmark_R(self):
        header = [
            'rate',
            'tool',
            'qc',
            'variable',
            'value',
        ]
        collection = []
        for record in self.collection:
            rate = record['simulated']
            for benchmark in record['benchmark']:
                tool = benchmark['tool']
                for qc in [ 'fail', 'pass', 'both' ]:
                # for qc in [ 'pass', 'fail', 'both' ]:
                    qnode = benchmark['multiplex']['barcode'][0]['noise'][qc]
                    for variable in [ 'TP', 'FP', 'FN', 'precision', 'recall', 'fscore' ]:
                        record = [
                            rate,
                            tool,
                            qc,
                            variable,
                            qnode[variable],
                        ]
                        collection.append(record)

        collection.sort(key=lambda i: i[2])
        collection.sort(key=lambda i: i[1])
        collection.sort(key=lambda i: i[0])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def summarize_decoder_with_noise_accuracy_benchmark_R(self):
        header = [
            'rate',
            'tool',
            'rank',
            'qc',
            'variable',
            'value',
        ]
        collection = []
        for record in self.collection:
            rate = record['simulated']
            for benchmark in record['benchmark']:
                tool = benchmark['tool']
                for rank in [ 'real' ]:
                # for rank in [ 'noise', 'real', 'both' ]:
                    for qc in [ 'pass', 'fail', 'both' ]:
                    # for qc in [ 'pass', 'fail', 'both' ]:
                        qnode = benchmark['multiplex']['decoder'][rank][qc]
                        for variable in [ 'FDR', 'MR', 'TP', 'FN', 'FP' ]:
                        # for variable in [ 'TP', 'FP', 'FN', 'FDR', 'MR', 'precision', 'recall' ]:
                            record = [
                                rate,
                                tool,
                                rank,
                                qc,
                                variable,
                                qnode[variable],
                            ]
                            collection.append(record)

        collection.sort(key=lambda i: i[3])
        collection.sort(key=lambda i: i[2])
        collection.sort(key=lambda i: i[1])
        collection.sort(key=lambda i: i[0])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def execute(self):
        self.collect_accuracy_benchmark()
        if self.instruction['preset'] == 'json':
            print(to_json(self.collection))

        elif self.instruction['preset'] == 'noise':
            self.summarize_noise_accuracy_benchmark()

        elif self.instruction['preset'] == 'multiplex':
            self.summarize_decoder_accuracy_benchmark()

        elif self.instruction['preset'] == 'multiplex_barcode':
            self.summarize_barcode_accuracy_benchmark()


        elif self.instruction['preset'] == 'quality_distribution_R':
            self.summarize_decoder_quality_distribution_R()

        elif self.instruction['preset'] == 'noise_R':
            self.summarize_noise_accuracy_benchmark_R()

        elif self.instruction['preset'] == 'multiplex_R':
            self.summarize_decoder_accuracy_benchmark_R()

        elif self.instruction['preset'] == 'multiplex_barcode_R':
            self.summarize_barcode_accuracy_benchmark_R()
