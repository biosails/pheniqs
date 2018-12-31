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

class Analyze(TranscodeSAM):
    def __init__(self, ontology):
        TranscodeSAM.__init__(self, ontology)

    def report(self):
        report = deepcopy(self.model)
        report['genealogy']['document type'] = 'analysis report'
        report['genealogy']['tool'] = self.instruction['tool']
        remove_compiled(report)
        if 'report' in self.instruction and self.instruction['report']:
            with io.open(self.instruction['report'], 'w') as file:
                file.write(to_json(report))

    def load_decoder_model(self, decoder):
        decoder['barcode by index.compiled'] = [ None ] * (decoder['barcode cardinality'] + 1)
        decoder['barcode by index.compiled'][0] = decoder['unclassified']
        for barcode in decoder['codec'].values():
            decoder['barcode by index.compiled'][barcode['index']] = barcode

        for qc in ['qc fail', 'qc pass']:
            for barcode in decoder['barcode by index.compiled']:
                barcode[qc] = {
                    'count': 0,
                    'FP': 0,
                    'FN': 0,
                    'TP': 0,
                    # 'TN': 0,
                    'precision': 0,
                    'recall': 0,
                    'FDR': 0,
                    'MR': 0,
                }

        if 'classified' not in decoder:
            decoder['classified'] = {}

        for qc in ['qc fail', 'qc pass']:
            decoder[qc] = {
                'count': 0,
                'FP': 0,
                'FN': 0,
                'TP': 0,
                # 'TN': 0,
                'precision': 0,
                'recall': 0,
                'FDR': 0,
                'MR': 0,
            }
            decoder['classified'][qc] = {
                'count': 0,
                'FP': 0,
                'FN': 0,
                'TP': 0,
                # 'TN': 0,
                'precision': 0,
                'recall': 0,
                'FDR': 0,
                'MR': 0,
            }

        decoder['macro'] = {}
        for qc in ['qc fail', 'qc pass']:
            decoder['macro'][qc] = {
                'precision': 0,
                'recall': 0,
                'FDR': 0,
                'MR': 0,
            }

        decoder['barcode by rg.compiled'] = {}
        decoder['unclassified']['RG'] = '{}:{}'.format(self.flowcell_id, 'undetermined')
        decoder['barcode by rg.compiled'][decoder['unclassified']['RG']] = decoder['unclassified']
        for barcode in decoder['codec'].values():
            barcode['RG'] = '{}:{}'.format(self.flowcell_id, ''.join(barcode['barcode']))
            decoder['barcode by rg.compiled'][barcode['RG']] = barcode

    def load_model(self):
        if TranscodeSAM.load_model(self):
            for topic in [ 'multiplex', 'cellular', 'molecular' ]:
                if topic in self.model:
                    if isinstance(self.model[topic], dict):
                        self.load_decoder_model(self.model[topic])
                        self.decoder_by_index.append(self.model[topic])

                    elif isinstance(self.model[topic], list):
                        for decoder in self.model[topic]:
                            self.load_decoder_model(decoder)
                            self.decoder_by_index.append(decoder)
                    return True
            else: raise NoConfigurationFileError('model file {} not found'.format(self.instruction['model']))
        return False

    def load(self):
        self.load_model()

    def is_qc_fail(self, read):
        flag = int(read['segment'][0]['fixed'][1])
        return (flag & 0x200 == 0x200)

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            qname = read['segment'][0]['fixed'][0]
            self.parse_qname(read)
            for hint, decoder in zip(read['hint'], self.decoder_by_index):

                true_barcode = decoder['barcode by index.compiled'][hint]
                if self.is_qc_fail(read):
                    qc = 'qc fail'
                else:
                    qc = 'qc pass'

                # if RG is not present assign to unclassified
                if 'RG' in read['segment'][0]['auxiliary']:
                    RG = read['segment'][0]['auxiliary']['RG']['VALUE']
                else:
                    RG = decoder['barcode by index.compiled'][0]['RG']

                decoded_barcode = decoder['barcode by rg.compiled'][RG]
                decoded_barcode[qc]['count'] += 1

                if true_barcode['RG'] == decoded_barcode['RG']:
                    true_barcode[qc]['TP'] += 1
                    # for barcode in self.barcode_by_index[decoder['key']]:
                    #     if barcode['index'] != true_barcode['index']:
                    #         barcode[qc]['TN'] += 1
                else:
                    true_barcode[qc]['FN'] += 1
                    decoded_barcode[qc]['FP'] += 1
                    # for barcode in self.barcode_by_index[decoder['key']]:
                    #     if barcode['index'] != true_barcode['index'] and barcode['index'] != decoded_barcode['index']:
                    #         barcode[qc]['TN'] += 1

    def flush(self):
        pass

    def finalize(self):
        for decoder in self.decoder_by_index:
            for barcode in decoder['barcode by index.compiled']:
                if barcode['index'] > 0:
                    for qc in ['qc fail', 'qc pass']:
                        if barcode[qc]['TP'] > 0 or barcode[qc]['FP'] > 0:
                            barcode[qc]['FDR'] = barcode[qc]['FP'] / (barcode[qc]['FP'] + barcode[qc]['TP'])
                            barcode[qc]['precision'] = barcode[qc]['TP'] / (barcode[qc]['FP'] + barcode[qc]['TP'])

                        if barcode[qc]['FN'] > 0 or barcode[qc]['TP'] > 0:
                            barcode[qc]['recall'] = barcode[qc]['TP'] / (barcode[qc]['TP'] + barcode[qc]['FN'])
                            barcode[qc]['MR'] = barcode[qc]['FN'] / (barcode[qc]['TP'] + barcode[qc]['FN'])

                        decoder[qc]['count'] += barcode[qc]['count']
                        decoder[qc]['FP'] += barcode[qc]['FP']
                        decoder[qc]['FN'] += barcode[qc]['FN']
                        decoder[qc]['TP'] += barcode[qc]['TP']
                        # decoder[qc]['TN'] += barcode[qc]['TN']

                        if barcode['index'] > 0:
                            decoder['classified'][qc]['count'] += barcode[qc]['count']
                            decoder['classified'][qc]['FP'] += barcode[qc]['FP']
                            decoder['classified'][qc]['FN'] += barcode[qc]['FN']
                            decoder['classified'][qc]['TP'] += barcode[qc]['TP']
                            # decoder['classified'][qc]['TN'] += barcode[qc]['TN']

                        decoder['macro'][qc]['FDR'] += barcode[qc]['FDR']
                        decoder['macro'][qc]['precision'] += barcode[qc]['precision']
                        decoder['macro'][qc]['recall'] += barcode[qc]['recall']
                        decoder['macro'][qc]['MR'] += barcode[qc]['MR']

            for qc in ['qc fail', 'qc pass']:
                if decoder[qc]['TP'] > 0 or decoder[qc]['FP'] > 0:
                    decoder[qc]['FDR'] = decoder[qc]['FP'] / (decoder[qc]['FP'] + decoder[qc]['TP'])
                    decoder[qc]['precision'] = decoder[qc]['TP'] / (decoder[qc]['FP'] + decoder[qc]['TP'])

                if decoder[qc]['FN'] > 0 or decoder[qc]['TP'] > 0:
                    decoder[qc]['recall'] = decoder[qc]['TP'] / (decoder[qc]['TP'] + decoder[qc]['FN'])
                    decoder[qc]['MR'] = decoder[qc]['FN'] / (decoder[qc]['TP'] + decoder[qc]['FN'])

                if decoder['classified'][qc]['TP'] > 0 or decoder['classified'][qc]['FP'] > 0:
                    decoder['classified'][qc]['FDR'] = decoder['classified'][qc]['FP'] / (decoder['classified'][qc]['FP'] + decoder['classified'][qc]['TP'])
                    decoder['classified'][qc]['precision'] = decoder['classified'][qc]['TP'] / (decoder['classified'][qc]['FP'] + decoder['classified'][qc]['TP'])

                if decoder['classified'][qc]['FN'] > 0 or decoder['classified'][qc]['TP'] > 0:
                    decoder['classified'][qc]['recall'] = decoder['classified'][qc]['TP'] / (decoder['classified'][qc]['TP'] + decoder['classified'][qc]['FN'])
                    decoder['classified'][qc]['MR'] = decoder['classified'][qc]['FN'] / (decoder['classified'][qc]['TP'] + decoder['classified'][qc]['FN'])

            for qc in ['qc fail', 'qc pass']:
                decoder['macro'][qc]['FDR'] /= decoder['barcode cardinality']
                decoder['macro'][qc]['precision'] /= decoder['barcode cardinality']
                decoder['macro'][qc]['recall'] /= decoder['barcode cardinality']
                decoder['macro'][qc]['MR'] /= decoder['barcode cardinality']
