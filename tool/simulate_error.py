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

class SimulateError(TranscodeSAM):
    def __init__(self, ontology):
        TranscodeSAM.__init__(self, ontology)

    @property
    def noise_model(self):
        return self.model['noise model']

    def load_model(self):
        if TranscodeSAM.load_model(self):
            self.model['genealogy']['error simulation id'] = str(uuid.uuid4())
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
        self.load_report()

    def report(self):
        if 'report' in self.instruction and self.instruction['report']:
            report = deepcopy(self.model)
            remove_compiled(report)
            report['genealogy']['document type'] = 'error simulation report';
            with io.open(self.instruction['report'], 'w') as file:
                file.write(to_json(report))

    def load_decoder_model(self, decoder):
        decoder['barcode by index.compiled'] = [ None ] * (decoder['barcode cardinality'] + 1)
        decoder['barcode by index.compiled'][0] = decoder['unclassified']
        for barcode in decoder['codec'].values():
            decoder['barcode by index.compiled'][barcode['index']] = barcode

        if 'transform' in decoder:
            decoder['transform.compiled'] = self.parse_rule(decoder['transform'])
            if 'substitution frequency transform' in decoder:
                decoder['substitution frequency transform.compiled'] = self.parse_rule(decoder['substitution frequency transform'])

        if 'error' in self.instruction and self.instruction['error']:
            decoder['calibrated'] = True
            decoder['original barcode error rate'] = decoder['expected barcode error']
            decoder['requested barcode error rate'] = self.instruction['error']
            normalization_factor = decoder['requested barcode error rate'] / decoder['original barcode error rate']
            self.log.info('calibrating expected error rate to %s with factor %s', decoder['requested barcode error rate'], normalization_factor)
            phred_probability_base = pow(10.0, -0.1)
            decoder['phred scale transform.compiled'] = []
            for q in range(64):
                o = pow(phred_probability_base, q)
                p = max(0, min(1, o * normalization_factor))
                v = round(-10.0 * math.log10(p))
                v = max(2, min(63, v))
                decoder['phred scale transform.compiled'].append(v)
        else:
            decoder['calibrated'] = False

        decoder['quality distribution'] = []
        for q in range(64):
            decoder['quality distribution'].append(0)

        self.load_substitution_model(decoder)

    def load_substitution_model(self, decoder):
        if 'substitution model' in decoder:
            decoder['substitution by quality.compiled'] = { 'segment': [] }
            for s in decoder['substitution model']['frequency']['segment']:
                segment = []
                for c in s['nucleotide substitution frequency']:
                    cycle = {}
                    for true, substitute in c.items():
                        nucleotide = []
                        for q in self.phred_scale:
                            correct = 1.0 - q['error']
                            quality = deepcopy(q)
                            quality['density'] = []

                            # no substitution
                            quality['density'].append({ 'nucleotide': true, 'density': correct })

                            # substitution
                            for o in substitute:
                                quality['density'].append({ 'nucleotide': o['to'], 'density': correct + o['density'] * q['error'] })

                            nucleotide.append(quality)
                        cycle[true] = nucleotide
                    segment.append(cycle)
                decoder['substitution by quality.compiled']['segment'].append(segment)
        else:
            raise BadConfigurationError('missing substitution model in decoder')

    def load_report(self):
        for decoder in self.decoder_by_index:
            decoder['expected event'] = 0
            decoder['expected barcode error compensation'] = 0
            decoder['expected barcode error'] = 0
            decoder['unclassified']['simulated substitution'] = 0
            decoder['simulated classified substitution'] = 0
            decoder['simulated substitution'] = 0
            for barcode in decoder['codec'].values():
                barcode['simulated substitution'] = 0

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            self.parse_qname(read)

            # convert segment SEQ and QUAL to a numpy char array
            for segment in read['segment']:
                segment['fixed'][9] = numpy.array(list(segment['fixed'][9]))
                segment['fixed'][10] = numpy.array(list(segment['fixed'][10]))

            for hint, decoder in zip(read['hint'], self.decoder_by_index):
                barcode = decoder['barcode by index.compiled'][hint]

                # simulate errors
                for token, barcode_segment_index in zip(decoder['transform.compiled']['token'], range(decoder['barcode segment cardinality'])):
                    for token_cycle, segment_cycle in zip(range(token['length']), range(token['start'], token['end'])):
                        nucleotide = read['segment'][token['segment']]['fixed'][9][segment_cycle]
                        quality = decode_phred(read['segment'][token['segment']]['fixed'][10][segment_cycle])

                        if decoder['calibrated']:
                            # calibrate the quality scores to match the expected error rate
                            quality = decoder['phred scale transform.compiled'][quality]
                            read['segment'][token['segment']]['fixed'][10][segment_cycle] = encode_phred(quality)

                        decoder['quality distribution'][quality] += 1

                        m = decoder['substitution by quality.compiled']['segment'][barcode_segment_index][token_cycle][nucleotide][quality]
                        y = m['error'] - decoder['expected barcode error compensation']
                        t = decoder['expected barcode error'] + y
                        compensation = (t - decoder['expected barcode error']) - y
                        decoder['expected barcode error'] = t

                        substitute = 'N'
                        event = numpy.random.random()
                        decoder['expected event'] += event
                        for option in m['density']:
                            if event < option['density']:
                                substitute = option['nucleotide']
                                break

                        if substitute != nucleotide:
                            barcode['simulated substitution'] += 1
                            read['segment'][token['segment']]['fixed'][9][segment_cycle] = substitute
                            self.log.debug('event %f : substitute %s with quality %d for %s at cycle %d of segment %d', event, nucleotide, quality, substitute, token_cycle, barcode_segment_index)

            # convert segment SEQ and QUAL back to string
            for segment in read['segment']:
                segment['fixed'][9] = ''.join(segment['fixed'][9])
                segment['fixed'][10] = ''.join(segment['fixed'][10])

            self.output_buffer.append(read)

    def finalize(self):
        for decoder in self.decoder_by_index:
            decoder['quality distribution count'] = 0
            decoder['quality distribution mean'] = 0
            decoder['quality distribution sd'] = 0
            for q in range(64):
                decoder['quality distribution count'] += decoder['quality distribution'][q]
                decoder['quality distribution mean'] += decoder['quality distribution'][q] * q
            decoder['quality distribution mean'] /= decoder['quality distribution count']

            for q in range(64):
                diff = q - decoder['quality distribution mean']
                decoder['quality distribution sd'] += decoder['quality distribution'][q] * diff * diff
            decoder['quality distribution sd'] /= decoder['quality distribution count']
            # del decoder['quality distribution']

            # unclassified
            unclassified = decoder['unclassified']
            unclassified['substitution rate per barcode'] = 0
            unclassified['substitution rate per nucleotide'] = 0

            if unclassified['count'] > 0:
                unclassified['substitution rate per barcode'] = unclassified['simulated substitution'] / unclassified['count']

            if unclassified['simulated nucleotide'] > 0:
                unclassified['substitution rate per nucleotide'] = unclassified['simulated substitution'] / unclassified['simulated nucleotide']

            # classified
            for barcode in decoder['codec'].values():
                barcode['substitution rate per barcode'] = 0
                barcode['substitution rate per nucleotide'] = 0
                decoder['simulated classified substitution'] += barcode['simulated substitution']

                if barcode['count'] > 0:
                    barcode['substitution rate per barcode'] = barcode['simulated substitution'] / barcode['count']

                if barcode['simulated nucleotide'] > 0:
                    barcode['substitution rate per nucleotide'] = barcode['simulated substitution'] / barcode['simulated nucleotide']

            # del decoder['expected barcode error compensation']
            decoder['simulated substitution'] = unclassified['simulated substitution'] + decoder['simulated classified substitution']

            decoder['substitution rate per nucleotide'] = 0
            if decoder['simulated nucleotide'] > 0:
                decoder['expected event'] /= decoder['simulated nucleotide']
                decoder['expected barcode error'] /= decoder['simulated nucleotide']
                decoder['substitution rate per nucleotide'] = decoder['simulated substitution'] / decoder['simulated nucleotide']

            decoder['substitution rate per classified nucleotide'] = 0
            if decoder['simulated classified nucleotide'] > 0:
                decoder['substitution rate per classified nucleotide'] = decoder['simulated classified substitution'] / decoder['simulated classified nucleotide']

            decoder['expected barcode error deviation'] = 1.0 - decoder['expected barcode error'] / decoder['substitution rate per nucleotide']
            self.log.info('substitution rate per nucleotide %s', decoder['substitution rate per nucleotide'])
            self.log.info('expected barcode error deviation %s', decoder['expected barcode error deviation'])
