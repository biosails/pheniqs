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
import random
import re
import logging
import os.path
import hashlib

PHIX_QUALITY_SAMPLE_PATH = 'HK5NHBGXX_phix_barcode_quality'
QUALITY_SAMPLE_PATH = 'HK5NHBGXX_barcode_quality'
PHIX_SAMPLE_PATH = 'HK5NHBGXX_phix'

configuration = {
    'model path': None,
    'read header pattern': '@NB501157:4:HK5NHBGXX:1:{}:{}:{} 1:N:0:{}\n',
    'phix read header pattern': '@NB501157:4:HK5NHBGXX:1:{}:{}:{}:PHIX 1:N:0:{}\n',
    'expression': {},
}
configuration['expression']['sam record'] = re.compile(
    r"""^
    (?P<QNAME>[!-?A-~]{1,254})\t
    (?P<FLAG>[0-9]+)\t
    (?P<RNAME>\*|[!-()+-<>-~][!-~]*)\t
    (?P<POS>[0-9]+)\t
    (?P<MAPQ>[0-9]+)\t
    (?P<CIGAR>\*|([0-9]+[MIDNSHPX=])+)\t
    (?P<RNEXT>\*|=|[!-()+-<>-~][!-~]*)\t
    (?P<PNEXT>[0-9]+)\t
    (?P<TLEN>[-+]?[0-9]+)\t
    (?P<SEQ>\*|[A-Za-z=.]+)\t
    (?P<QUAL>[!-~]*)
    (?P<OPT>(\t([A-Za-z][A-Za-z0-9]:[AifZHB]:[^\s]+))+)?
    $""",
    re.VERBOSE)
configuration['expression']['sam optional field'] = re.compile(
    '(?P<TAG>[A-Za-z][A-Za-z0-9]):(?P<TYPE>[AifZHB]):(?P<VALUE>[!-~]+)')
configuration['expression']['sam optional field value'] = {
    'A': re.compile('^[!-~]$'),
    'i': re.compile('^[-+]?[0-9]+$'),
    'f': re.compile('^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'),
    'Z': re.compile('^[ !-~]+$'),
    'H': re.compile('^[0-9A-F]+$'),
    'B': re.compile('^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+$'),
}
configuration['expression']['phix record'] = re.compile(r'.*PHIX$')

# simulated phix reads. 
# wgsim 0.3.2
# wgsim -e 0.01 -N 100000 -1 8 -2 8 phix_NC_001422.fa phix1.fastq phix2.fastq

def to_json(node):
    print(json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4))

def load_report(log):
    if os.path.isfile(configuration['report path']):
        with io.open(configuration['report path'], 'r') as handle:
            report = json.loads(handle.read())
    else:
        report = { 'record': {}, }
    return report

def ammend_report(log, node):
    dirty = False
    report = load_report(log)
    for decoder in ['mdd','pamld']:
        if decoder in node:
            record = node[decoder]
            if record['sha1'] not in report['record']:
                report['record'][record['sha1']] = record
                dirty = True

    if dirty:
        with io.open(configuration['report path'], 'wb') as handle:
           handle.write(json.dumps(report, sort_keys=True, ensure_ascii=False, indent=4).encode('utf8'))

def csv(log):
    content = io.open('/dev/stdin', 'r').read()
    report = json.loads(content)
    head = [
        'decoder',
        'confidence',
        'variable',
        'value',
        'nature'
    ]
    print('\t'.join(head))

    for key, record in report['record'].items():
        if record['confidence'] is None:
            record['confidence'] = 0.0

        for variable in ['precision', 'recall', 'f score']:
        # for variable in ['accuracy', 'precision', 'recall', 'f score']:
            print('\t'.join([
                record['decoder'],
                '{0:.8f}'.format(1.0 - record['confidence']),
                variable.replace(' ', '_'),
                str(record[variable]),
                record['nature'],
            ]))

def csv_no_undetermined(log):
    content = io.open('/dev/stdin', 'r').read()
    report = json.loads(content)
    buffer = []
    for key, record in report['record'].items():
        record['FN'] = 0
        record['FP'] = 0
        record['TN'] = 0
        record['TP'] = 0

        if record['confidence'] is None:
            record['confidence'] = 0.99

        for barcode, library in record['barcode'].items():
            if '=' not in barcode:
                record['FN'] += library['FN']
                record['FP'] += library['FP']
                record['TN'] += library['TN']
                record['TP'] += library['TP']

        record['FDR'] = record['FP'] / (record['FP'] + record['TP'])
        record['precision'] = record['TP'] / (record['FP'] + record['TP'])
        record['recall'] = record['TP'] / (record['TP'] + record['FN'])
        record['MR'] = record['FN'] / (record['TP'] + record['FN'])

    for key, record in report['record'].items():
        buffer.append([
            record['decoder'].upper(),
            '{0:.8f}'.format(1.0 - record['confidence']),
            'precision',
            str(record['FDR']),
            record['nature'],
        ])
        buffer.append([
            record['decoder'].upper(),
            '{0:.8f}'.format(1.0 - record['confidence']),
            'recall',
            str(record['MR']),
            record['nature'],
        ])

    head = [
        'decoder',
        'confidence',
        'variable',
        'value',
        'nature'
    ]
    print('\t'.join(head))

    buffer.sort(key=lambda x: x[1])
    buffer.sort(key=lambda x: x[0])
    for row in buffer:
        print('\t'.join(row))

def fdr(log):
    content = io.open('/dev/stdin', 'r').read()
    report = json.loads(content)
    buffer = []
    for key, record in report['record'].items():
        record['FN'] = 0
        record['FP'] = 0
        record['TN'] = 0
        record['TP'] = 0

        for barcode, library in record['barcode'].items():
            if '=' not in barcode:
                record['FN'] += library['FN']
                record['FP'] += library['FP']
                record['TN'] += library['TN']
                record['TP'] += library['TP']

        record['FDR'] = record['FP'] / (record['FP'] + record['TP'])
        record['precision'] = record['TP'] / (record['FP'] + record['TP'])
        record['recall'] = record['TP'] / (record['TP'] + record['FN'])
        record['MR'] = record['FN'] / (record['TP'] + record['FN'])

    mdd = None
    for key, record in report['record'].items():
        if 'decoder' in record and record['decoder'] == 'mdd':
            mdd = record

    for key, record in report['record'].items():
        if 'decoder' in record and record['decoder'] == 'pamld':
            if(record['MR'] > mdd['MR']):
                record['DMR'] = -(100 * ((record['MR'] / mdd['MR']) - 1.0))
            else:
                record['DMR'] = 100 * ((mdd['MR'] / record['MR']) - 1.0)
            record['DMR'] = '{0:.2f}'.format(record['DMR'])

            if(record['FDR'] > 0):
                if(record['FDR'] > mdd['FDR']):
                    record['DFDR'] = -(100 * ((record['FDR'] / mdd['FDR']) - 1.0))
                else:
                    record['DFDR'] = 100 * ((mdd['FDR'] / record['FDR']) - 1.0)
                record['DFDR'] = '{0:.2f}'.format(record['DFDR'])

            else:
                record['DFDR'] = 'INF'

            row = [
                '{0:.8f}'.format(record['confidence']),
                str(record['TP']),
                str(record['FP']),
                str(record['FN']),
                '{0:.8f}'.format(round(record['FDR'], 8)),
                '{0:.8f}'.format(round(record['MR'], 8)),
                record['DFDR'],
                record['DMR'],
            ]
            buffer.append(row)

    head = [
        'confidence',
        'TP',
        'FP',
        'FN',
        'FDR',
        'MR',
        'Delta FDR',
        'Delta DMR',
    ]
    print('\t& '.join(head))

    buffer.sort(key=lambda x: x[0])
    for row in buffer:
        print('\t& '.join(row))

    row = [
        'NA',
        str(mdd['TP']),
        str(mdd['FP']),
        str(mdd['FN']),
        '{0:.8f}'.format(round(mdd['FDR'], 8)),
        '{0:.8f}'.format(round(mdd['MR'], 8)),
        'NA',
        'NA',
    ]
    print('\t& '.join(row))


def csv_only_undetermined(log):
    content = io.open('/dev/stdin', 'r').read()
    report = json.loads(content)
    head = [
        'decoder',
        'confidence',
        'variable',
        'value',
        'nature'
    ]
    print('\t'.join(head))

    for key, record in report['record'].items():
        precision = 0
        recall = 0
        fscore = 0
        count = 0
        if record['confidence'] is None:
            record['confidence'] = 0.99

        for barcode, library in record['barcode'].items():
            if '=' in barcode:
                precision = library['precision']
                recall = library['recall']

        fscore = f_score(precision, recall)
        confidence = 1.0 - record['confidence']

        print('\t'.join([
            record['decoder'],
            '{0:.8f}'.format(confidence),
            'precision',
            str(precision),
            record['nature'],
        ]))

        print('\t'.join([
            record['decoder'],
            '{0:.8f}'.format(confidence),
            'recall',
            str(recall),
            record['nature'],
        ]))

        print('\t'.join([
            record['decoder'],
            '{0:.8f}'.format(confidence),
            'f_score',
            str(fscore),
            record['nature'],
        ]))

def load_model(log):
    content = io.open(configuration['model path']).read()
    model = json.loads(content)
    for quality in model['quality'].values():
        for mer in quality['count'].values():
            mer['order'] = sorted(mer['count'].keys())
            mer['size'] = len(mer['order'])
            mer['accumulative'] = []
            total = 0
            for index,key in enumerate(mer['order']):
                total += mer['count'][key]
                mer['accumulative'].append(total)

    model['barcode length'] = len(model['barcode'][0])
    model['barcode size'] = len(model['barcode'])
    model['null barcode'] = '=' * model['barcode length']
    model['phix concentration'] = 1.0 / model['phix frequency']
    model['library concentration'] = (1.0  - model['phix concentration']) / float(model['barcode size'])

    model['tile'] = model['tile start']
    model['x'] = model['x start']
    model['y'] = model['y start']
    configuration['model'] = model

def load_sample(log):
    sample = []
    content = io.open(configuration['model']['quality sample path']).readlines()
    for line in content:
        line = line.strip()
        if line:
            record = line.split(' ')
            record[0] = int(record[0])
            for i in range(record[0]):
                sample.append(record[1])

    configuration['sample size'] = len(sample)
    configuration['sample'] = sample

def load_phix_sample(log):
    if 'phix quality sample path' in configuration['model']:
        sample = []
        content = io.open(configuration['model']['phix quality sample path']).readlines()
        for line in content:
            line = line.strip()
            if line:
                record = line.split(' ')
                record[0] = int(record[0])
                for i in range(record[0]):
                    sample.append(record[1])

        configuration['phix sample size'] = len(sample)
        configuration['phix sample'] = sample

def load_phix(log):
    sample = []
    content = io.open(configuration['model']['phix sample path']).readlines()
    for line in content:
        line = line.strip()
        if line:
            sample.append(line)
    configuration['phix size'] = len(sample)
    configuration['phix'] = sample

def pick_barcode(length=None):
    position = random.randrange(configuration['model']['barcode size'])
    value = configuration['model']['barcode'][position]
    if length is not None:
        value = value[:length]
    return value

def pick_phix(length=None):
    position = random.randrange(configuration['phix size'])
    value = configuration['phix'][position]
    if length is not None:
        value = value[:length]
    return value

def pick_sample(length=None):
    position = random.randrange(configuration['sample size'])
    value = configuration['sample'][position]
    if length is not None:
        value = value[:length]
    return value

def pick_phix_sample(length=None):
    position = random.randrange(configuration['phix sample size'])
    value = configuration['phix sample'][position]
    if length is not None:
        value = value[:length]
    return value

def pick_output(quality, truth):
    if quality in configuration['model']['quality']:
        q = configuration['model']['quality'][quality]
        if truth in q['count']:
            t = q['count'][truth]
            event = random.randrange(t['total'])

            position = 0
            while(event > t['accumulative'][position]):
                position += 1
            return t['order'][position]
    else:
        return truth

def simulate_record(log, handle):
    if configuration['model']['tile'] < configuration['model']['tile end']:
        configuration['model']['tile'] += 1
    else:
        configuration['model']['tile'] = configuration['model']['tile start']

        if configuration['model']['y'] < configuration['model']['y end']:
            configuration['model']['y'] += 1
        else:
            configuration['model']['y'] = configuration['model']['y start']
            configuration['x'] += 1

    phix = random.randrange(configuration['model']['phix frequency'])
    length = configuration['model']['barcode length']
    if phix > 0:
        sample = pick_sample(length)
        barcode = pick_barcode(length)
        output = []
        check = []
        for base in zip(sample, barcode):
            o = pick_output(base[0], base[1])
            output.append(o)
            if o == base[1]:
                check.append('-')
            else:
                check.append('*')
        check = ''.join(check)
        output = ''.join(output)
        handle.write(configuration['read header pattern'].format (
            configuration['model']['tile'],
            configuration['model']['x'],
            configuration['model']['y'], 
            '').encode('ascii'))
    else:
        sample = pick_phix_sample(length)
        barcode = pick_phix(length)
        output = []
        check = []
        for base in zip(sample, barcode):
            o = pick_output(base[0], base[1])
            output.append(o)
            if o == base[1]:
                check.append('-')
            else:
                check.append('*')
        check = ''.join(check)
        output = ''.join(output)
        handle.write(configuration['phix read header pattern'].format (
            configuration['model']['tile'], 
            configuration['model']['x'], 
            configuration['model']['y'],
            '').encode('ascii'))

    handle.write('{}{}\n'.format(output, barcode).encode('ascii'))
    handle.write('+\n'.encode('ascii'))
    handle.write('{}{}\n'.format(sample, sample).encode('ascii'))

def simulate(log):
    load_model(log)
    load_sample(log)
    load_phix(log)
    load_phix_sample(log)
    left = configuration['model']['population size']
    with io.open('/dev/stdout', 'wb') as handle:
        while(left > 0):
            simulate_record(log, handle)
            left -= 1

def fill(log, capacity):
    collection = []
    for line in sys.stdin:
        if line:
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
            #
            # 2  RNAME   string     Reference sequence NAME
            # 3  POS     int        1-based leftmost mapping POSition
            # 4  MAPQ    int        MAPping Quality
            # 5  CIGAR   string     CIGAR string
            # 6  RNEXT   string     Reference name of the mate/next read
            # 7  PNEXT   int        Position of the mate/next read
            # 8  TLEN    int        observed Template LENgth
            # 9  SEQ     string     segment SEQuence
            # 10 QUAL    string     Phred QUALity+33
            #
            # MD    Z   String for mismatching positions. [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)
            # AS    i   Alignment score generated by aligner
            # NM    i   Edit distance to the reference, including ambiguous bases but excluding clipping
            # XS        Suboptimal alignment score
            match = configuration['expression']['sam record'].search(line)
            if match:
                sam_record = match.groupdict()

                # Parsing
                for field in [
                    'FLAG',
                    'POS',
                    'MAPQ',
                    'PNEXT',
                    'TLEN',
                ]:
                    sam_record[field] = int(sam_record[field])

                if 'OPT' in sam_record:
                    optional = sam_record['OPT'].strip('\t').split('\t')
                    for o in optional:
                        m = configuration['expression']['sam optional field'].search(o)
                        if m:
                            o = m.groupdict()
                            if configuration['expression']['sam optional field value'][o['TYPE']].match(o['VALUE']):
                                if o['TYPE'] is 'i':
                                    o['VALUE'] = int(o['VALUE'])
                                elif o['TYPE'] is 'f':
                                    o['VALUE'] = float(o['VALUE'])
                                sam_record[o['TAG']] = o['VALUE']
                            else:
                                log.error('ignoring invalid %s optional field %s', o['TYPE'], o['VALUE'])

                if configuration['expression']['phix record'].match(sam_record['QNAME']):
                    sam_record['PHIX'] = True
                else:
                    sam_record['PHIX'] = False

                for field in [
                    'OPT',
                    'CIGAR',
                    'FLAG',
                    'MAPQ',
                    'PNEXT',
                    'POS',
                    'RNAME',
                    'RNEXT',
                    'TLEN',
                    'QUAL',
                    'RG'
                ]:
                    del sam_record[field]

                if 'DQ' not in sam_record:
                    sam_record['DQ'] = 0
                collection.append(sam_record)
            else:
                log.error('invalid sam syntax %s', line)
        if len(collection) >= capacity:
            break
    return collection

def f_score(precision, recall):
    if(precision + recall) > 0:
        return 2.0 * (precision * recall) / (precision + recall)
    else:
        return 0

def checksum(decoder, nature, confidence=None):
    checksum = hashlib.sha1()
    checksum.update(decoder.encode('utf8'))
    checksum.update(nature.encode('utf8'))
    if confidence is not None:
        checksum.update('{0:.10f}'.format(round(confidence, 10)).encode('utf8'))
    return checksum.hexdigest()

def analyze(log):
    load_model(log)
    model = configuration['model']
    node = {
        'count': 0,
        'barcode': list(model['barcode']),
        'nature': model['nature'],
        'confidence': configuration['confidence'],
        'pamld': {
            'sha1': checksum('pamld', model['nature'], configuration['confidence']),
            'nature': model['nature'],
            'confidence': configuration['confidence'],
            'decoder': 'pamld',
            'wrong': 0,
            'correct': 0,
            'accuracy': 0.0,
            'precision': 0.0,
            'recall': 0.0,
            'f score': 0.0,
            'barcode': {},
        },
        'mdd': {
            'sha1': checksum('mdd', model['nature']),
            'nature': model['nature'],
            'confidence': None,
            'decoder': 'mdd',
            'wrong': 0,
            'correct': 0,
            'accuracy': 0.0,
            'precision': 0.0,
            'recall': 0.0,
            'f score': 0.0,
            'barcode': {},
        },
    }
    node['barcode'].append(model['null barcode'])
    for decoder in ['mdd','pamld']:
        measure = node[decoder]
        for barcode in node['barcode']:
            measure['barcode'][barcode] = {
                'TP': 0,
                'TN': 0,
                'FP': 0,
                'FN': 0,
            }

    collection = fill(log, 65536)
    while(collection):
        for sam_record in collection:
            analyze_record(log, node, sam_record)
        collection = fill(log, 65536)

    if node['count'] > 0:
        for decoder in ['mdd','pamld']:
            measure = node[decoder]
            measure['accuracy'] = float(measure['correct']) / float(node['count'])
            for key, barcode in measure['barcode'].items():
                value = barcode['TP'] + barcode['FP']
                if value:
                    barcode['precision'] = barcode['TP'] / value
                else:
                    barcode['precision'] = 0.0

                value = barcode['TP'] + barcode['FN']
                if value:
                    barcode['recall'] = barcode['TP'] / value
                else:
                    barcode['recall'] = 0.0
                barcode['f score'] = f_score(barcode['precision'], barcode['recall'])

                if '=' in key:
                    measure['precision'] += (model['phix concentration'] * barcode['precision'])
                    measure['recall'] += (model['phix concentration'] * barcode['recall'])
                else:
                    measure['precision'] += (model['library concentration'] * barcode['precision'])
                    measure['recall'] += (model['library concentration'] * barcode['recall'])
            measure['f score'] = f_score(measure['precision'], measure['recall'])

    ammend_report(log, node)
    to_json(node)

def analyze_record(log, node, sam_record):
    ID = sam_record['QNAME']        # read id
    CC = sam_record['SEQ']          # correct barcode
    PL = sam_record['XL']           # PAMLD barcode
    MD = sam_record['XM']           # MDD barcode
    PX = sam_record['PHIX']         # is the read a simulated phix
    # BC = sam_record['BC']           # observed barcode
    # DQ = float(sam_record['DQ'])    # PAMLD decoding error probability
    # EE = float(sam_record['EE'])    # Expected number of errors in barcode

    if PX: CC = configuration['model']['null barcode']
    if CC in node['barcode']:
        node['count'] += 1
        for decoder in ['mdd', 'pamld']:
            measure = node[decoder]

            if decoder == 'mdd': DB = MD
            else: DB = PL

            if CC in measure['barcode']:
                if DB == CC:            
                    measure['correct'] += 1
                    measure['barcode'][CC]['TP'] += 1
                    for k,v in measure['barcode'].items():
                        if k != CC: v['TN'] += 1
                else:
                    measure['wrong'] += 1
                    measure['barcode'][CC]['FN'] += 1
                    measure['barcode'][DB]['FP'] += 1
                    for k,v in measure['barcode'].items():
                        if k != CC and k != DB: v['TN'] += 1
    else:
        log.warning('unknown barcode {} in read {}'.format(CC, ID))

def learn(log):
    head = False
    model = { 
        'quality sample path': QUALITY_SAMPLE_PATH,
        'phix sample path': PHIX_SAMPLE_PATH,
        'phix quality sample path': PHIX_QUALITY_SAMPLE_PATH,
        'phix frequency': 100,
        'population size': 10000000,
        'quality': {},
        'count': 0,
        'barcode': [],
        'tile start': 1143,
        'tile end': 16347,
        'x start': 11101,
        'x end': 15000,
        'y start': 10000,
        'y end': 15000,
    }
    for line in sys.stdin:
        if line:
            if not head:
                head = True
            else:
                row = line.strip().split()
                DX = int(row[1])
                OBS = row[8]
                OBQ = row[9]
                PBS = row[11]
                DBS = row[13]
                DBD = int(row[14])

                # enriched model, only reads that have at least one error
                # if '=' not in PBS and '=' not in DBS and DX == 0 and DBD > 0:
                #
                # natural model, all reads
                if '=' not in PBS and '=' not in DBS and DX == 0:

                    BC = DBS
                    if BC not in model['barcode']:
                        model['barcode'].append(BC)

                    model['count'] += 1

                    # print(DX, OBS, OBQ, PBS, DBS, DBD)
                    for cycle, nucleotide in enumerate(BC):
                        score = OBQ[cycle]
                        converted = OBS[cycle]
                        
                        if score not in model['quality']:
                            model['quality'][score] = {
                                'score': ord(score) - 33,
                                'total': 0,
                                'count': {},
                            }
                        quality = model['quality'][score]
                        quality['total'] += 1 

                        if nucleotide not in quality['count']:
                            quality['count'][nucleotide] = {
                                'total': 0,
                                'count': {},
                            }
                        conversion = quality['count'][nucleotide]
                        conversion['total'] += 1 
                        if converted not in conversion['count']:
                            conversion['count'][converted] = 0

                        conversion['count'][converted] += 1

    for quality in model['quality'].values():
        for nucleotide in quality['count'].values():
            nucleotide['conversion'] = {}
            for k,v in nucleotide['count'].items():
                nucleotide['conversion'][k] = v / nucleotide['total']
    to_json(model)

def print_help():
    print('simulate.py <verb> <model path>\nverb: learn, analyze, simulate, csv')

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    log = logging.getLogger('simulate')

    if len(sys.argv) > 1:
        action = sys.argv[1]
        if action == 'fdr':
            fdr(log)

        elif action == 'csv':
            # csv(log)
            csv_no_undetermined(log)
            # csv_only_undetermined(log)

        elif len(sys.argv) > 2:
            configuration['model path'] = sys.argv[2]

            if action == 'learn':
                learn(log)

            elif action == 'analyze':
                if len(sys.argv) > 4:
                    configuration['report path'] = sys.argv[3]
                    configuration['confidence'] = float(sys.argv[4])
                    analyze(log)
                else:
                    print('simulate.py analyze <model path> <report path> <confidence>')

            elif action == 'simulate':
                simulate(log)

            else:
                print('unknown action {}'.format(action))
                print_help()
        else:
            print_help()
    else:
        print_help()

if __name__ == '__main__':
    main()

