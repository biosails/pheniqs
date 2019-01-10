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

from core.error import *
from core import Job
from core import to_json

class Summarize(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)

    @property
    def experiment(self):
        return self.ontology['experiment']

    def collect_speed_benchmark(self):
        header = [
            'requested',
            'expected',
            'simulated',
            'req/sim',
            'exp/sim',
            'tool',
            'kernel',
            'user',
            'real',
            'memory',
            'cpu',
        ]
        collection = []
        # collection.append(header)
        for substitution in self.experiment['substitution'].values():
            if 'model' in substitution:
                requested = substitution['model']['multiplex']['requested barcode error rate']
                expected = substitution['model']['multiplex']['expected barcode error']
                simulated = substitution['model']['multiplex']['substitution rate per nucleotide']

                if 'deml demultiplex execution' in substitution:
                    execution = substitution['deml demultiplex execution']
                    cpu = (100 * (execution['kernel'] + execution['user']) / execution['real'])
                    record = [
                        'deml',
                        execution['kernel'],
                        execution['user'],
                        execution['real'],
                        execution['memory'],
                        cpu,
                        requested,
                        # expected,
                        simulated,
                        requested / simulated,
                        # expected / simulated,
                    ]
                    collection.append(record)

                if 'pheniqs demultiplex execution' in substitution:
                    execution = substitution['pheniqs demultiplex execution']
                    cpu = (100 * (execution['kernel'] + execution['user']) / execution['real'])
                    record = [
                        'pheniqs',
                        execution['kernel'],
                        execution['user'],
                        execution['real'],
                        execution['memory'],
                        cpu,
                        requested,
                        # expected,
                        simulated,
                        requested / simulated,
                        # expected / simulated,
                    ]
                    collection.append(record)

                if 'sense prior execution' in substitution:
                    execution = substitution['sense prior execution']
                    cpu = (100 * (execution['kernel'] + execution['user']) / execution['real'])
                    record = [
                        'prior',
                        execution['kernel'],
                        execution['user'],
                        execution['real'],
                        execution['memory'],
                        cpu,
                        requested,
                        # expected,
                        simulated,
                        requested / simulated,
                        # expected / simulated,
                    ]
                    collection.append(record)

        collection.sort(key=lambda i: i[0])
        collection.sort(key=lambda i: i[6])
        content = '\n'.join([','.join([str(field) for field in record]) for record in collection])
        print(content)

    def collect_classified_accuracy(self):
        def summarize_barcode(decoder, barcode):
            for qc in ['qc fail', 'qc pass']:
                for item in [ 'count', 'TP', 'FP', 'FN' ]:
                    decoder[qc][item] += barcode[qc][item]

        def summarize_decoder(decoder):
            for qc in ['qc fail', 'qc pass']:
                decoder[qc]['count'] = 0

            for barcode in decoder['codec'].values():
                summarize_barcode(decoder, barcode)

            for qc in ['qc fail', 'qc pass']:
                TP = decoder[qc]['TP']
                FP = decoder[qc]['FP']
                FN = decoder[qc]['FN']
                if TP > 0 or FP > 0:
                    decoder[qc]['precision']    = TP / (FP + TP)
                    decoder[qc]['FDR']          = FP / (FP + TP)
                if FN > 0 or TP > 0:
                    decoder[qc]['recall']       = TP / (TP + FN)
                    decoder[qc]['MR']           = FN / (TP + FN)

        header = [
            'tool',
            'simulated',
            'qc',
            'count',
            'TP',
            'FP',
            'FN',
            'FDR',
            'MR',
            'precision',
            'recall',
        ]
        collection = []

        # collection.append(header)
        for substitution in self.experiment['substitution'].values():
            if 'model' in substitution:
                requested = substitution['model']['multiplex']['requested barcode error rate']
                expected = substitution['model']['multiplex']['expected barcode error']
                simulated = substitution['model']['multiplex']['substitution rate per nucleotide']

                for tool in [ 'deml', 'pheniqs' ]:
                    akey = '{} demultiplex analysis'.format(tool)
                    if akey in substitution:
                        analysis = substitution[akey]
                        summarize_decoder(analysis['multiplex'])
                        for qc in [ 'pass', 'fail' ]:
                            qnode = analysis['multiplex']['qc {}'.format(qc)]
                            record = [
                                tool,
                                simulated,
                                qc,
                                qnode['count'],
                                qnode['TP'],
                                qnode['FP'],
                                qnode['FN'],
                                qnode['FDR'],
                                qnode['MR'],
                                qnode['precision'],
                                qnode['recall'],
                            ]
                            if qc == 'pass': collection.append(record)
                            # collection.append(record)

        collection.sort(key=lambda i: i[0])
        collection.sort(key=lambda i: i[2])
        collection.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def collect_accuracy_benchmark(self):
        def summarize_barcode(decoder, barcode):
            for qc in ['qc fail', 'qc pass']:
                for item in [ 'count', 'TP', 'FP', 'FN' ]:
                    decoder[qc][item] += barcode[qc][item]
                    if barcode['index'] > 0:
                        decoder['classified'][qc][item] += barcode[qc][item]

        def summarize_decoder(decoder):
            for qc in ['qc fail', 'qc pass']:
                decoder[qc]['count'] = 0

            summarize_barcode(decoder, decoder['unclassified'])
            for barcode in decoder['codec'].values():
                summarize_barcode(decoder, barcode)

            for qc in ['qc fail', 'qc pass']:
                TP = decoder[qc]['TP']
                FP = decoder[qc]['FP']
                FN = decoder[qc]['FN']
                if TP > 0 or FP > 0:
                    decoder[qc]['precision']    = TP / (FP + TP)
                    decoder[qc]['FDR']          = FP / (FP + TP)
                if FN > 0 or TP > 0:
                    decoder[qc]['recall']       = TP / (TP + FN)
                    decoder[qc]['MR']           = FN / (TP + FN)

                TP = decoder['classified'][qc]['TP']
                FP = decoder['classified'][qc]['FP']
                FN = decoder['classified'][qc]['FN']
                if TP > 0 or FP > 0:
                    decoder['classified'][qc]['precision']  = TP / (FP + TP)
                    decoder['classified'][qc]['FDR']        = FP / (FP + TP)

                if FN > 0 or TP > 0:
                    decoder['classified'][qc]['recall']     = TP / (TP + FN)
                    decoder['classified'][qc]['MR']         = FN / (TP + FN)

        header = [
            'tool',
            'simulated',
            'qc',
            'count',
            'TP',
            'FP',
            'FN',
            'FDR',
            'MR',
            'precision',
            'recall',
        ]
        collection = []

        # collection.append(header)
        for substitution in self.experiment['substitution'].values():
            if 'model' in substitution:
                requested = substitution['model']['multiplex']['requested barcode error rate']
                expected = substitution['model']['multiplex']['expected barcode error']
                simulated = substitution['model']['multiplex']['substitution rate per nucleotide']

                for tool in [ 'deml', 'pheniqs' ]:
                    akey = '{} demultiplex analysis'.format(tool)
                    if akey in substitution:
                        analysis = substitution[akey]
                        summarize_decoder(analysis['multiplex'])
                        qsum = [
                            tool,
                            simulated,
                            'all',
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                        ]
                        for qc in [ 'pass', 'fail' ]:
                            qnode = analysis['multiplex']['qc {}'.format(qc)]
                            qsum[3] += qnode['count']
                            qsum[4] += qnode['TP']
                            qsum[5] += qnode['FP']
                            qsum[6] += qnode['FN']

                            record = [
                                tool,
                                simulated,
                                qc,
                                qnode['count'],
                                qnode['TP'],
                                qnode['FP'],
                                qnode['FN'],
                                qnode['FDR'],
                                qnode['MR'],
                                qnode['precision'],
                                qnode['recall'],
                            ]
                            if qc == 'pass': collection.append(record)
                            # collection.append(record)

                        # FDR =         FP / ( FP + TP )
                        # Precision =   TP / ( FP + TP )
                        # MR =          FN / ( FN + TP )
                        # Recall =      TP / ( FN + TP )
                        TP = qsum[4]
                        FP = qsum[5]
                        FN = qsum[6]
                        if TP > 0 or FP > 0:
                            qsum[7]  = FP / (FP + TP) # FDR
                            qsum[9]  = TP / (FP + TP) # Precision

                        if FN > 0 or TP > 0:
                            qsum[10] = TP / (FN + TP) # Recall
                            qsum[8]  = FN / (FN + TP) # MR
                        # collection.append(qsum)

        collection.sort(key=lambda i: i[0])
        collection.sort(key=lambda i: i[2])
        collection.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in record]) for record in collection]))

    def summarize_accuracy_benchmark(self):
        def summarize_barcode(n, barcode):
            d = n['decoder']
            b = n['barcode'][barcode['index']]
            for qc in ['fail', 'pass']:
                for item in [ 'count', 'TP', 'FP', 'FN' ]:
                    value = barcode['qc {}'.format(qc)][item]
                    b[qc][item] = value
                    b['both'][item] += value

                    # update decoder
                    d['all'][qc][item] += value
                    d['all']['both'][item] += value
                    if barcode['index'] > 0:
                        d['classified'][qc][item] += value
                        d['classified']['both'][item] += value

            for qc in ['fail', 'pass', 'both']:
                if b[qc]['TP'] > 0 or b[qc]['FP'] > 0:
                    b[qc]['precision']    = b[qc]['TP'] / (b[qc]['FP'] + b[qc]['TP'])
                    b[qc]['FDR']          = b[qc]['FP'] / (b[qc]['FP'] + b[qc]['TP'])
                if b[qc]['FN'] > 0 or b[qc]['TP'] > 0:
                    b[qc]['recall']       = b[qc]['TP'] / (b[qc]['TP'] + b[qc]['FN'])
                    b[qc]['MR']           = b[qc]['FN'] / (b[qc]['TP'] + b[qc]['FN'])

        def summarize_decoder(n, decoder):
            n['decoder'] = {}
            n['barcode'] = []

            d = n['decoder']
            for rank in [ 'all', 'classified' ]:
                d[rank] = {}
                for qc in ['fail', 'pass', 'both']:
                    d[rank][qc] = {}
                    for item in [ 'count', 'TP', 'FP', 'FN' ]:
                        d[rank][qc][item] = 0

            for index in range(decoder['barcode cardinality'] + 1):
                b = { 'index': index }
                for qc in ['fail', 'pass', 'both']:
                    b[qc] = {}
                    for item in [ 'count', 'TP', 'FP', 'FN' ]:
                        b[qc][item] = 0
                n['barcode'].append(b)

            summarize_barcode(n, decoder['unclassified'])
            for barcode in decoder['codec'].values():
                summarize_barcode(n, barcode)

            for rank in [ 'all', 'classified' ]:
                for qc in ['fail', 'pass', 'both']:
                    for item in [ 'count', 'TP', 'FP', 'FN' ]:
                        if d[rank][qc]['TP'] > 0 or d[rank][qc]['FP'] > 0:
                            d[rank][qc]['precision']    = d[rank][qc]['TP'] / (d[rank][qc]['FP'] + d[rank][qc]['TP'])
                            d[rank][qc]['FDR']          = d[rank][qc]['FP'] / (d[rank][qc]['FP'] + d[rank][qc]['TP'])
                        if d[rank][qc]['FN'] > 0 or d[rank][qc]['TP'] > 0:
                            d[rank][qc]['recall']       = d[rank][qc]['TP'] / (d[rank][qc]['TP'] + d[rank][qc]['FN'])
                            d[rank][qc]['MR']           = d[rank][qc]['FN'] / (d[rank][qc]['TP'] + d[rank][qc]['FN'])

        collection = []
        for substitution in self.experiment['substitution'].values():
            if 'model' in substitution:
                record = {
                    'requested': substitution['model']['multiplex']['requested barcode error rate'],
                    'expected': substitution['model']['multiplex']['expected barcode error'],
                    'simulated': substitution['model']['multiplex']['substitution rate per nucleotide'],
                    'benchmark': []
                }
                for tool in [ 'deml', 'pheniqs' ]:
                    akey = '{} demultiplex analysis'.format(tool)
                    if akey in substitution:
                        analysis = { 'tool': tool, }
                        if 'multiplex' in substitution[akey]:
                            analysis['multiplex'] = {}
                            summarize_decoder(analysis['multiplex'], substitution[akey]['multiplex']),
                        record['benchmark'].append(analysis)
                collection.append(record)
        print(to_json(collection))



    def collect_read_count(self):
        def summarize_barcode(decoder, barcode):
            for qc in ['qc fail', 'qc pass']:
                decoder[qc]['count'] += barcode[qc]['count']

        def summarize_decoder(decoder):
            for qc in ['qc fail', 'qc pass']:
                decoder[qc]['count'] = 0

            summarize_barcode(decoder, decoder['unclassified'])
            for barcode in decoder['codec'].values():
                summarize_barcode(decoder, barcode)

        for substitution in self.experiment['substitution'].values():
            if 'model' in substitution:
                simulated = substitution['model']['multiplex']['substitution rate per nucleotide']
                for tool in [ 'deml', 'pheniqs' ]:
                    akey = '{} demultiplex analysis'.format(tool)
                    if akey in substitution:
                        analysis = substitution[akey]
                        summarize_decoder(analysis['multiplex'])
                        record = [
                            simulated,
                            tool,
                            analysis['multiplex']['qc pass']['count'],
                            analysis['multiplex']['qc fail']['count'],
                            analysis['multiplex']['qc fail']['count'] + analysis['multiplex']['qc pass']['count']
                        ]
                        print(','.join([str(field) for field in record]))

    def execute(self):
        # self.collect_speed_benchmark()
        # self.collect_accuracy_benchmark()
        # self.collect_classified_accuracy()
        # self.collect_read_count()
        self.summarize_accuracy_benchmark()

    # def make_summary(self):
    #     self.ontology['summary'] = { 'total': {}, 'barcode': {} }
    #     for record in self.ontology['experiment']:
    #         rate = record['rate']
    #         tool = record['tool']
    #         self.log.info('collecting data from %s', record['basename'])
    #         with io.open(record['path'], 'rb') as file:
    #             experiment = json.loads(file.read().decode('utf8'))
    #             if rate not in self.summary['total']: self.summary['total'][rate] = {}
    #             if tool not in self.summary['total'][rate]: self.summary['total'][rate][tool] = {}
    #             for qc in  self.scale['qc']:
    #                 self.summary['total'][rate][tool][qc] = experiment['decoder']['A5KVK'][qc]
    #
    #             for k,v in experiment['decoder']['A5KVK']['codec'].items():
    #                 if k not in self.summary['barcode']: self.summary['barcode'][k] = {}
    #                 if rate not in self.summary['barcode'][k]: self.summary['barcode'][k][rate] = {}
    #                 if tool not in self.summary['barcode'][k][rate]: self.summary['barcode'][k][rate][tool] = {}
    #                 for qc in  self.scale['qc']:
    #                     self.summary['barcode'][k][rate][tool][qc] = v[qc]
    #
    #     for rate in self.summary['total'].values():
    #         for tool in rate.values():
    #             both = {}
    #             for qc in tool.values():
    #                 for k,v in qc.items():
    #                     if k not in both: both[k] = v
    #                     else: both[k] += v
    #             tool['both'] = both
    #
    #     for barcode in self.summary['barcode'].values():
    #         for rate in barcode.values():
    #             for tool in rate.values():
    #                 both = {}
    #                 for qc in tool.values():
    #                     for k,v in qc.items():
    #                         if k not in both: both[k] = v
    #                         else: both[k] += v
    #                 tool['both'] = both
    #
    # def write_summary_csv_for_node(self, name, node):
    #     output_path = '{}_summary.csv'.format(name)
    #     content = []
    #     head = [ 'rate' ]
    #     for qc in self.header['qc'] + ['B']:
    #         for field in self.header['field']:
    #             for tool in self.header[ 'tool' ]:
    #                 head.append('{}:{}:{}'.format(tool, qc, field))
    #     content.append(','.join(head))
    #
    #     for rate in sorted(node.keys()):
    #         row = [ '{}'.format(rate) ]
    #         for qc in self.scale['qc'] + ['both']:
    #             for field in self.scale['field']:
    #                 for tool in self.scale['tool']:
    #                     row.append(node[rate][tool][qc][field])
    #         if row:
    #             content.append(','.join([str(i) for i in row]))
    #
    #     with io.open(output_path, 'w') as file:
    #         file.write('\n'.join(content))
    #
    # def write_csv_for_each_barcode(self):
    #     self.write_summary_csv_for_node('total', self.summary['total'])
    #     for barcode in sorted(self.summary['barcode'].keys()):
    #         self.write_summary_csv_for_node(barcode,self.summary['barcode'][barcode])
    #
    # def write_summary_csv(self):
    #     output_path = 'summary.csv'
    #     content = []
    #     head = [ 'barcode', 'rate' ]
    #     for qc in self.header['qc']:
    #         for field in self.header['field']:
    #             for tool in self.header[ 'tool' ]:
    #                 head.append('{}:{}:{}'.format(tool, qc, field))
    #     content.append(','.join(head))
    #
    #     for barcode in sorted(self.summary['barcode'].keys()):
    #         for rate in sorted(self.summary['barcode'][barcode].keys()):
    #             row = [ barcode, '{}'.format(rate) ]
    #             for qc in self.scale['qc']:
    #                 for field in self.scale['field']:
    #                     for tool in self.scale['tool']:
    #                         row.append(self.summary['barcode'][barcode][rate][tool][qc][field])
    #             content.append(','.join([str(i) for i in row]))
    #
    #     with io.open(output_path, 'w') as file:
    #         file.write('\n'.join(content))
    #
    # def execute(self):
    #     self.load()
    #
    #     self.write_csv_for_each_barcode()
    #     self.write_summary_csv()
    #
    #     self.finalize()
    #
    # def collect_per_barcode():
    #     pass
    #
    # def collect_decoder_analysis(self, decoder):
    #     meta_decoder = {}
    #     for item in [
    #         'key',
    #         'index',
    #         'noise',
    #         'macro',
    #         'qc fail',
    #         'qc pass',
    #         'classified count',
    #         'simulated classified nucleotide',
    #         'simulated classified substitution',
    #         'simulated noise',
    #         'simulated noise deviation',
    #         'simulated nucleotide',
    #         'simulated substitution',
    #         'substitution rate per classified nucleotide',
    #         'substitution rate per nucleotide',
    #     ]:
    #         if item in decoder:
    #             meta_decoder[item] = decoder[item]
    #
    #     if 'codec' in decoder:
    #         meta_decoder['codec'] = {}
    #         for barcode in decoder['codec'].values():
    #             meta_barcode = {}
    #             meta_decoder['codec'][barcode['key']] = meta_barcode
    #             for item in [
    #                 'key',
    #                 'RG'
    #                 'index',
    #                 'barcode',
    #                 'concentration',
    #                 'count',
    #                 'qc fail',
    #                 'qc pass',
    #                 'simulated concentration',
    #                 'simulated nucleotide',
    #                 'simulated substitution',
    #                 'substitution rate per barcode',
    #                 'substitution rate per nucleotide'
    #             ]:
    #                 if item in barcode:
    #                     meta_barcode[item] = barcode[item]
    #
    #     if 'unclassified' in decoder:
    #         meta_decoder['unclassified'] = {}
    #         for item in [
    #             'RG'
    #             'index',
    #             'concentration',
    #             'count',
    #             'qc fail',
    #             'qc pass',
    #             'simulated concentration',
    #             'simulated nucleotide',
    #             'simulated substitution',
    #             'substitution rate per barcode',
    #             'substitution rate per nucleotide'
    #         ]:
    #             if item in decoder['unclassified']:
    #                 meta_decoder['unclassified'][item] = decoder['unclassified'][item]
    #     return meta_decoder
    #
    # def execute(self):
    #     self.load()
    #     self.collect_analysis()
    #     self.finalize()
    #     self.persist()
    #
    # def finalize(self):
    #     pass
