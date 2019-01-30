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
import logging
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
        self.log = logging.getLogger('Analyze')
        self.genealogy['tool'] = self.instruction['tool']
        self.model

        for term in [
            'noise model',
        ]:
            if term in self.model:
                del self.model[term]

        if 'multiplex' in self.model:
            for term in [
                'substitution frequency transform',
                'substitution model',
            ]:
                if term in self.model['multiplex']:
                    del self.model['multiplex'][term]

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

        self.instruction['input'] = os.path.join(self.home, self.instruction['input'])

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

                # Noise accumulator:
                #   TP : a noise read was correctly classified as noise
                #   FN :
                #       noise: a noise read was incorrectly classified to a real barcode
                #   FP :
                #       noise: a real read was classified as noise
                #
                # Real barcode accumulator:
                #   TP : real read is correctly classified
                #   FN :
                #       real: a read from this barcode was incorrectly classified to another real barcode
                #       noise: a read from this barcode was incorrectly classified as noise
                #   FP :
                #       real: a real read from another barcode was incorrectly classified to this barcode
                #       noise: a noise read was incorrectly classified to this barcode
                if true_barcode['index'] > 0:
                    # real read
                    if true_barcode['index'] == decoded_barcode['index']:
                        # real read is correctly classified
                        true_barcode['accumulate']['real'][qc]['TP'] += 1               # TP real
                    else:
                        # real read is incorrectly classified
                        if decoded_barcode['index'] > 0:
                            # real read is incorrectly classified to another real barcode
                            # true barcode is real
                            # decoded barcode is real
                            true_barcode['accumulate']['real'][qc]['FN'] += 1           # FN real
                            decoded_barcode['accumulate']['real'][qc]['FP'] += 1        # FP real
                        else:
                            # real read is classified as noise
                            # true barcode is real
                            # decoded barcode is noise
                            true_barcode['accumulate']['noise'][qc]['FN'] += 1          # FN real
                            decoded_barcode['accumulate']['noise'][qc]['FP'] += 1       # FP noise
                else:
                    # noise read
                    if decoded_barcode['index'] > 0:
                        # noise read is incorrectly classified to a real barcode
                        # true barcode is noise
                        # decoded barcode is real
                        true_barcode['accumulate']['noise'][qc]['FN'] += 1              # FN noise
                        decoded_barcode['accumulate']['noise'][qc]['FP'] += 1           # FP real
                    else:
                        # noise read is correctly classified as noise
                        true_barcode['accumulate']['noise'][qc]['TP'] += 1              # TP noise

    def flush(self):
        pass

    def finalize(self):
        pass
