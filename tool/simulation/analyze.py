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

import os
from core.error import *
from simulation.transcode import Transcode

# FDR = FP / ( FP + TP )
# Precision = TP / ( FP + TP )
# MR = FN / ( FN + TP )
# recall = TP / ( FN + TP )

def is_qc_fail(read):
    flag = int(read['segment'][0]['fixed'][1])
    return (flag & 0x200 == 0x200)

class Analyze(Transcode):
    def __init__(self, ontology):
        Transcode.__init__(self, ontology)
        self.genealogy['tool'] = self.instruction['tool']

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    def load(self):
        def load_barcode_model(barcode):
            barcode['accumulate'] = {}
            for rank in ['real', 'noise']:
                barcode['accumulate'][rank] = {}
                for qc in ['fail', 'pass']:
                    barcode['accumulate'][rank][qc] = {}
                    for item in [ 'count', 'TP', 'FP', 'FN' ]:
                        barcode['accumulate'][rank][qc][item] = 0

        def load_decoder_model(decoder):
            # load enumerated barcode by index
            decoder['barcode by index.compiled'] = [ None ] * (decoder['barcode cardinality'] + 1)
            decoder['barcode by index.compiled'][0] = decoder['unclassified']
            for barcode in decoder['codec'].values():
                decoder['barcode by index.compiled'][barcode['index']] = barcode

            for barcode in decoder['barcode by index.compiled']:
                load_barcode_model(barcode)

        def load_multiplex_model(decoder):
            # load RG to barcode mapping
            decoder['barcode by rg.compiled'] = {}

            decoder['unclassified']['RG'] = '{}:{}'.format(self.flowcell_id, 'undetermined')
            decoder['barcode by rg.compiled'][decoder['unclassified']['RG']] = decoder['unclassified']

            for barcode in decoder['codec'].values():
                barcode['RG'] = '{}:{}'.format(self.flowcell_id, ''.join(barcode['barcode']))
                decoder['barcode by rg.compiled'][barcode['RG']] = barcode

        Transcode.load(self)

        for topic in [ 'multiplex', 'cellular', 'molecular' ]:
            if topic in self.model:
                if isinstance(self.model[topic], dict):
                    load_decoder_model(self.model[topic])
                    self.decoder_by_index.append(self.model[topic])

                elif isinstance(self.model[topic], list):
                    for decoder in self.model[topic]:
                        load_decoder_model(decoder)
                        self.decoder_by_index.append(decoder)

        if 'multiplex' in self.model:
            load_multiplex_model(self.model['multiplex'])

        self.instruction['input'] = os.path.join(self.home, self.location['{} demultiplex path'.format(self.genealogy['tool'])])

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            qname = read['segment'][0]['fixed'][0]
            self.parse_qname(read)

            if 'multiplex' in self.model:
                hint = read['hint'][self.model['multiplex']['index']]

                # find the true barcode from the header
                true_barcode = self.model['multiplex']['barcode by index.compiled'][hint]

                # find the decoded barcode from the assigned RG
                if 'RG' in read['segment'][0]['auxiliary']:
                    decoded_barcode = self.model['multiplex']['barcode by rg.compiled'][read['segment'][0]['auxiliary']['RG']['VALUE']]
                else:
                    # if RG is not present assign to unclassified
                    decoded_barcode = self.model['multiplex']['barcode by rg.compiled'][self.model['multiplex']['barcode by index.compiled'][0]['RG']]

                # fail or pass quality control
                if is_qc_fail(read):
                    qc = 'fail'
                else:
                    qc = 'pass'

                if true_barcode['index'] > 0:
                    # real read, not noise
                    if true_barcode['index'] == decoded_barcode['index']:
                        # read is correctly classified
                        true_barcode['accumulate']['real'][qc]['TP'] += 1
                    else:
                        # read is incorrectly classified
                        if decoded_barcode['index'] > 0:
                            # read is classified to another, real, barcode
                            true_barcode['accumulate']['real'][qc]['FN'] += 1
                            decoded_barcode['accumulate']['real'][qc]['FP'] += 1
                        else:
                            # read is classified as noise
                            true_barcode['accumulate']['noise'][qc]['FN'] += 1
                            decoded_barcode['accumulate']['noise'][qc]['FP'] += 1
                else:
                    # noise read : true_barcode['index'] == 0
                    if decoded_barcode['index'] > 0:
                        # noise read was classified to a real barcode
                        decoded_barcode['accumulate']['noise'][qc]['FP'] += 1
                        true_barcode['accumulate']['noise'][qc]['FN'] += 1
                    else:
                        # noise read is correctly classified as noise
                        true_barcode['accumulate']['noise'][qc]['TP'] += 1

    def flush(self):
        pass

    def finalize(self):
        pass
        # for decoder in self.decoder_by_index:
        #     for barcode in decoder['barcode by index.compiled']:
        #         for qc in ['qc fail', 'qc pass']:
        #             TP = barcode[qc]['TP']
        #             FP = barcode[qc]['FP']
        #             FN = barcode[qc]['FN']
        #             if TP > 0 or FP > 0:
        #                 barcode[qc]['precision']    = TP / (FP + TP)
        #                 barcode[qc]['FDR']          = FP / (FP + TP)
        #             if TP > 0 or FN > 0:
        #                 barcode[qc]['recall']       = TP / (TP + FN)
        #                 barcode[qc]['MR']           = FN / (TP + FN)
        #
        #             decoder[qc]['count'] += barcode[qc]['count']
        #             decoder[qc]['FP'] += FP
        #             decoder[qc]['FN'] += FN
        #             decoder[qc]['TP'] += TP
        #             if barcode['index'] > 0:
        #                 decoder['classified'][qc]['count'] += barcode[qc]['count']
        #                 decoder['classified'][qc]['FP'] += FP
        #                 decoder['classified'][qc]['FN'] += FN
        #                 decoder['classified'][qc]['TP'] += TP
        #             decoder['macro'][qc]['precision'] += barcode[qc]['precision']
        #             decoder['macro'][qc]['recall'] += barcode[qc]['recall']
        #             decoder['macro'][qc]['FDR'] += barcode[qc]['FDR']
        #             decoder['macro'][qc]['MR'] += barcode[qc]['MR']
        #
        #     for qc in ['qc fail', 'qc pass']:
        #         TP = decoder[qc]['TP']
        #         FP = decoder[qc]['FP']
        #         FN = decoder[qc]['FN']
        #         if TP > 0 or FP > 0:
        #             decoder[qc]['precision']    = TP / (FP + TP)
        #             decoder[qc]['FDR']          = FP / (FP + TP)
        #         if FN > 0 or TP > 0:
        #             decoder[qc]['recall']       = TP / (TP + FN)
        #             decoder[qc]['MR']           = FN / (TP + FN)
        #
        #         TP = decoder['classified'][qc]['TP']
        #         FP = decoder['classified'][qc]['FP']
        #         FN = decoder['classified'][qc]['FN']
        #         if TP > 0 or FP > 0:
        #             decoder['classified'][qc]['precision']  = TP / (FP + TP)
        #             decoder['classified'][qc]['FDR']        = FP / (FP + TP)
        #
        #         if FN > 0 or TP > 0:
        #             decoder['classified'][qc]['recall']     = TP / (TP + FN)
        #             decoder['classified'][qc]['MR']         = FN / (TP + FN)
