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

from sam import *

def print_json(node):
    def handler(o):
        if isinstance(o, numpy.ndarray):
            return ''.join(o)

        elif isinstance(o, Decimal):
            return '{:10.17f}'.format(o)

        return None

    return print(json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4, default=handler))

def f_score(precision, recall):
    if(precision + recall) > 0:
        return 2.0 * (precision * recall) / (precision + recall)
    else:
        return 0

class Analyze(SamTranscode):
    def __init__(self):
        SamTranscode.__init__(self, 'analyze')
        self.ontology['barcode by index'] = {}
        self.ontology['barcode by RG'] = {}

    @property
    def barcode_by_index(self):
        return self.ontology['barcode by index']

    @property
    def barcode_by_RG(self):
        return self.ontology['barcode by RG']

    def load_RG(self):
        for decoder in self.decoder_by_index:
            self.barcode_by_RG[decoder['key']] = {}
            decoder['unclassified']['RG'] = '{}:{}'.format(self.flowcell_id, 'undetermined')
            self.barcode_by_RG[decoder['key']][decoder['unclassified']['RG']] = decoder['unclassified']

            for barcode in decoder['codec'].values():
                barcode['RG'] = '{}:{}'.format(self.flowcell_id, ''.join(barcode['barcode']))
                self.barcode_by_RG[decoder['key']][barcode['RG']] = barcode

    def load(self):
        if os.path.exists(self.instruction['path']):
            self.instruction['working directory'] = os.path.dirname(self.instruction['path'])
            self.log.debug('loading simulation %s', self.instruction['path'])
            with io.open(self.instruction['path'], 'rb') as file:
                self.ontology['preset'] = json.loads(file.read().decode('utf8'))

            self.decoder_by_index = [ None ] * len(self.preset['decoder'])
            for key,decoder in self.preset['decoder'].items():
                # self.token_model[key] = { 'token': self.parse_token(decoder['transform']['token']) }
                self.barcode_by_index[key] = [ None ] * (decoder['barcode cardinality'] + 1)
                self.barcode_by_index[key][0] = decoder['unclassified']
                for barcode in decoder['codec'].values():
                    self.barcode_by_index[key][barcode['index']] = barcode

                for barcode in self.barcode_by_index[key]:
                    barcode['FP'] = 0
                    barcode['FN'] = 0
                    barcode['TP'] = 0
                    barcode['TN'] = 0
                    barcode['FDR'] = 0
                    barcode['precision'] = 0
                    barcode['recall'] = 0
                    barcode['MR'] = 0

                decoder['correct'] = 0
                decoder['wrong'] = 0
                decoder['macro'] = {
                    'FDR': 0,
                    'precision': 0,
                    'recall': 0,
                    'MR': 0,
                }
                decoder['FP'] = 0
                decoder['FN'] = 0
                decoder['TP'] = 0
                decoder['TN'] = 0
                decoder['FDR'] = 0
                decoder['precision'] = 0
                decoder['recall'] = 0
                decoder['MR'] = 0
                self.decoder_by_index[decoder['index']] = decoder

                self.log.info('loaded %s decoder', key)
            self.load_RG()
        else: raise NoConfigurationFileError('simulation {} not found'.format(self.instruction['path']))

    def parse_qname(self, read):
        # PHENIQS:0:A5KVK:0000016164:S49:C
        # 0:1:2:3:4:5
        # 0 instrument id
        # 1 run count
        # 2 flowcell id
        # 3 read counter
        # 4 - separated true sample barcode indexs
        # 5 - separated true cellular barcode indexs
        read['fact'] = {}
        parsed = read['segment'][0]['fixed'][0].split(':')
        if len(parsed) > 4 and len(parsed[4]) > 1:
            read['fact']['sample'] = [int(i) for i in parsed[4].strip('S').split('-')]

        if len(parsed) > 5 and len(parsed[5]) > 1:
            read['fact']['cellular'] = [int(i) for i in parsed[5].strip('C').split('-')]

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            self.parse_qname(read)
            if 'sample' in read['fact']:
                if 'RG' in read['segment'][0]['auxiliary']:
                    RG = read['segment'][0]['auxiliary']['RG']['VALUE']
                    for decoder in self.decoder_by_index:
                        true_barcode_index = read['fact']['sample'][decoder['index']]
                        true_barcode = self.barcode_by_index[decoder['key']][true_barcode_index]
                        decoded_barcode = self.barcode_by_RG[decoder['key']][RG]
                        if true_barcode['RG'] == decoded_barcode['RG']:
                            decoder['correct'] += 1
                            true_barcode['TP'] += 1
                            for barcode in self.barcode_by_index[decoder['key']]:
                                if barcode['index'] != true_barcode['index']:
                                    barcode['TN'] += 1
                        else:
                            decoder['wrong'] += 1
                            true_barcode['FN'] += 1
                            decoded_barcode['FP'] += 1
                            for barcode in self.barcode_by_index[decoder['key']]:
                                if barcode['index'] != true_barcode['index'] and barcode['index'] != decoded_barcode['index']:
                                    barcode['TN'] += 1
            else:
                self.log.warning('read %s is missing a read group', read[segment][0]['fixed'][0])

            # self.output_buffer.append(read)

    def flush(self):
        pass

    def finalize(self):
        for decoder in self.decoder_by_index:
            for barcode in self.barcode_by_index[decoder['key']]:
                if barcode['TP'] > 0 or barcode['FP'] > 0:
                    barcode['FDR'] = barcode['FP'] / (barcode['FP'] + barcode['TP'])
                    barcode['precision'] = barcode['TP'] / (barcode['FP'] + barcode['TP'])

                if barcode['FN'] > 0 or barcode['TP'] > 0:
                    barcode['recall'] = barcode['TP'] / (barcode['TP'] + barcode['FN'])
                    barcode['MR'] = barcode['FN'] / (barcode['TP'] + barcode['FN'])

                decoder['FP'] += barcode['FP']
                decoder['FN'] += barcode['FN']
                decoder['TP'] += barcode['TP']
                decoder['TN'] += barcode['TN']
                decoder['macro']['FDR'] += barcode['FDR']
                decoder['macro']['precision'] += barcode['precision']
                decoder['macro']['recall'] += barcode['recall']
                decoder['macro']['MR'] += barcode['MR']

            if decoder['TP'] > 0 or decoder['FP'] > 0:
                decoder['FDR'] = decoder['FP'] / (decoder['FP'] + decoder['TP'])
                decoder['precision'] = decoder['TP'] / (decoder['FP'] + decoder['TP'])

            if decoder['FN'] > 0 or decoder['TP'] > 0:
                decoder['recall'] = decoder['TP'] / (decoder['TP'] + decoder['FN'])
                decoder['MR'] = decoder['FN'] / (decoder['TP'] + decoder['FN'])

            decoder['macro']['FDR'] /= (decoder['barcode cardinality'] + 1)
            decoder['macro']['precision'] /= (decoder['barcode cardinality'] + 1)
            decoder['macro']['recall'] /= (decoder['barcode cardinality'] + 1)
            decoder['macro']['MR'] /= (decoder['barcode cardinality'] + 1)

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)

    pipeline = None

    try:
        pipeline = Analyze()
        pipeline.execute()

    except (
        PermissionDeniedError,
        NoOverwriteError,
        DownloadError,
        CommandFailedError,
        NoConfigurationFileError,
        BadConfigurationError,
        UnsupportedError,
        SequenceError
    ) as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except(KeyboardInterrupt, SystemExit) as e:
        if e.code != 0:
            logging.getLogger('main').critical(e)
            sys.exit(1)

    finally:
        if pipeline: pipeline.close()

    sys.exit(0)

if __name__ == '__main__':
    main()
