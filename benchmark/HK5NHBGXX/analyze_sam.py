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

# Kraken version 0.10.6-unreleased
# Copyright 2013-2015, Derrick Wood (dwood@cs.jhu.edu)

# dustmasker: 1.0.0
# Package: blast 2.2.31, build Jun  2 2015 10:20:04

# Program: bwa (alignment via Burrows-Wheeler transformation)
# Version: 0.7.8-r455

import re
import io
import json
import numpy
import sys
import logging
from io import StringIO, BytesIO
from subprocess import Popen, PIPE

# default flowcell
FLOWCELL = 'HK5NHBGXX'

# number of processors on the machine
THREADS = 16

# number of reads in each batch.
# will stay in RAM during the processing so don't overdo this
BUFFER_SIZE = 8192 * 32

# paths to executables
DUSTMASKER = 'dustmasker'
KRAKEN = 'kraken'
KRAKEN_FILTER = 'kraken-filter'
BWA = 'bwa'

# path to directory with indexed references
REFERENCE = '/Users/lg/reference'

# path to taxonomy json file
TAXONOMY = '/Users/lg/code/pheniqs/benchmark/taxonomy.json'

# path to kraken db
KRAKEN_DB = '/volume/albireo/waverly/kraken/standard'

expression = {}
expression['sam record'] = re.compile(
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
    (?P<QUAL>[!-~]+)
    (?P<OPT>(\t([A-Za-z][A-Za-z0-9]:[AifZHB]:[^\s]+))+)?
    $""",
    re.VERBOSE
)
expression['sam optional field'] = re.compile('(?P<TAG>[A-Za-z][A-Za-z0-9]):(?P<TYPE>[AifZHB]):(?P<VALUE>[^:]+)')
expression['sam optional field value'] = {
    'A': re.compile('^[!-~]$'),
    'i': re.compile('^[-+]?[0-9]+$'),
    'f': re.compile('^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'),
    'Z': re.compile('^[ !-~]+$'),
    'H': re.compile('^[0-9A-F]+$'),
    'B': re.compile('^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+$'),
}
configuration = {
    'flowcell': {
        'HK5NHBGXX' : {
            'capacity': BUFFER_SIZE,
            'barcode': {
                '================': { 'LN': [0,        ], 'OG': [ 'PX'                   ]},
                'AAGAGGCAAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AAGAGGCAAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AAGAGGCAATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AAGAGGCACTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AAGAGGCATACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AAGAGGCATATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AAGAGGCATCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AAGAGGCATCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AGGCAGAAAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AGGCAGAAAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AGGCAGAAATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AGGCAGAACTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AGGCAGAATACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AGGCAGAATATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AGGCAGAATCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'AGGCAGAATCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CAGAGAGGAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CAGAGAGGAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CAGAGAGGATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CAGAGAGGCTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CAGAGAGGTACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CAGAGAGGTATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CAGAGAGGTCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CAGAGAGGTCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGAGGCTGAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGAGGCTGAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGAGGCTGATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGAGGCTGCTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGAGGCTGTACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGAGGCTGTATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGAGGCTGTCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGAGGCTGTCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGTACTAGAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGTACTAGAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGTACTAGATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGTACTAGCTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGTACTAGTACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGTACTAGTATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGTACTAGTCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CGTACTAGTCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CTCTCTACAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CTCTCTACAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CTCTCTACATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CTCTCTACCTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CTCTCTACTACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CTCTCTACTATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CTCTCTACTCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'CTCTCTACTCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GCTACGCTAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GCTACGCTAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GCTACGCTATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GCTACGCTCTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GCTACGCTTACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GCTACGCTTATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GCTACGCTTCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GCTACGCTTCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GGACTCCTAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GGACTCCTAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GGACTCCTATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GGACTCCTCTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GGACTCCTTACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GGACTCCTTATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GGACTCCTTCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GGACTCCTTCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GTAGAGGAAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GTAGAGGAAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GTAGAGGAATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GTAGAGGACTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GTAGAGGATACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GTAGAGGATATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GTAGAGGATCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'GTAGAGGATCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAAGGCGAAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAAGGCGAAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAAGGCGAATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAAGGCGACTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAAGGCGATACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAAGGCGATATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAAGGCGATCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAAGGCGATCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAGGCATGAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAGGCATGAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAGGCATGATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAGGCATGCTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAGGCATGTACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAGGCATGTATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAGGCATGTCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TAGGCATGTCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TCCTGAGCAGAGGATA': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TCCTGAGCAGGCTTAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TCCTGAGCATAGAGAG': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TCCTGAGCCTCCTTAC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TCCTGAGCTACTCCTT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TCCTGAGCTATGCAGT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TCCTGAGCTCTACTCT': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
                'TCCTGAGCTCTTACGC': { 'LN': [1, 2, 3, 4], 'OG': [ 'SC', 'CO', 'SP', 'NC' ]},
            },
            'organism': [
                'SC',
                'CO',
                'SP',
                'NC',
                'PX',
            ]
        }
    },
    'organism': {
        'SC': {
            'name': 'Saccharomyces cerevisiae',
            'path': REFERENCE +'/S288C_reference_sequence_R64-2-1_20150113.fsa'
        },
        'CO': {
            'name': 'Candida orthopsilosis',
            'path': REFERENCE + '/GCF_000315875.1_ASM31587v1_genomic.fna'
        },
        'SP': {
            'name': 'Saccharomyces cerevisiae YJM1527 plasmid 2 micron',
            'path': REFERENCE + '/CP004565.fa'
        },
        'NC': {
            'name': 'Naumovozyma castellii',
            'path': REFERENCE + '/GCF_000237345.1_ASM23734v1_genomic.fna'
        },
        'PX': {
            'name': 'phiX174 sensu lato',
            'path': REFERENCE +'/phix_NC_001422.fa'
        },
    }
}
def to_json(node):
    print(json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4))

def hamming(one, two):
    distance = 0
    for o, t in zip(one, two):
      if o != t:
          distance += 1
    return distance

class Pipeline(object):
    def __init__(self, configuration, flowcell):
        self.log = logging.getLogger('Pipeline')
        self.configuration = configuration
        self.flowcell = flowcell
        self.collection = None
        self.height = len(self.barcode)
        self.width = max([len(b) for b in self.barcode.keys()])
        self.barcodeSequence = sorted(self.barcode.keys())
        self.lowComplexity = re.compile('[atcgn]')

        # hamming distance and likelihood
        self.shortest_distance = self.width
        self.hammigDistance = numpy.zeros((self.height, self.height), dtype=numpy.int32)
        for i in range(self.height):
            for j in range(self.height):
                if i != j:
                    d = hamming(self.barcodeSequence[i],self.barcodeSequence[j])
                    self.hammigDistance[i][j] = d
                    self.shortest_distance = min(self.shortest_distance, d)

        self.randomLikelihood = 1 / pow(4, self.width)
        self.randomBarcodeLikelihood = len(self.barcode) / pow(4, self.width)
        self.randomLikelihoodCorrectableBarcode = self.height / pow(4, self.width - self.shortest_distance + 1)

        for barcode in self.barcode.values():
            barcode['LN'] = set(barcode['LN'])
            barcode['OG'] = set(barcode['OG'])

        self.division = {
            0: { 'key': 'BC', 'name': 'bacteria' },
            3: { 'key': 'PH', 'name': 'phages' },
            9: { 'key': 'VR', 'name': 'viruses' },
        }
        self.taxonomy = {}
        node = json.loads(io.open(TAXONOMY).read())
        for division_id, division in node.items():
            d_id = int(division_id)
            for taxonomy_id in division:
                self.taxonomy[taxonomy_id] = d_id

    @property
    def size(self):
        return len(self.collection)

    @property
    def full(self):
        return self.size >= self.capacity

    @property
    def empty(self):
        return self.size == 0

    @property
    def barcode(self):
        return self.flowcell['barcode']

    @property
    def organism(self):
        return self.flowcell['organism']

    @property
    def capacity(self):
        return self.flowcell['capacity']

    def run(self):
        self.head()
        self.cycle()

    def cycle(self):
        while self.fill():
            self.dust()

            # this will set the OG tag on bacterial and viral sequences 
            self.kraken()
            # self.filtered_kraken()

            # this will set the OG tag on the associated organisms
            for o in self.organism:
                self.bwa(o)
            self.analyze()
            self.csv()

    def json(self):
        to_json(self.collection)

    def head(self):
        print(' '.join([
            'CC',
            'DX',
            'WL',
            'WO',
            'OG',
            'AS',
            'MQ',
            'QS',
            'OBS',
            'OBQ',
            'BEE',
            'PBS',
            'PBD',
            'DBS',
            'DBD',
            'P',
            'CP',
            'REE',
            'RS',
            'RQ',
            'ID',
            'BS',
            'BD',
        ]))

    def csv(self):
        for r in self.collection.values():
            row = [
                r['CC'],    # 0
                r['DX'],    # 1
                r['WL'],    # 2
                r['WO'],    # 3
                r['OG'],    # 4
                r['AS'],    # 5
                r['MQ'],    # 6
                r['QS'],    # 7
                r['OBS'],   # 8
                r['OBQ'],   # 9
                r['BEE'],   # 10
                r['PBS'],   # 11
                r['PBD'],   # 12
                r['DBS'],   # 13
                r['DBD'],   # 14
                r['P'],     # 15
                r['CP'],    # 16
                r['REE'],   # 17
                r['RS'],    # 18
                r['RQ'],    # 19
                r['ID'],    # 20
                r['BS'],    # 21
                r['BD'],    # 22
            ]
            print(' '.join([str(c) for c in row]))

    def analyze(self):
        for record in self.collection.values():
            b = self.barcode[record['PBS']]
            if record['OG'] != b['OG']:
                record['WO'] = 1

            if record['WL'] > 0 or record['WO'] > 0:
                record['CC'] = 0

    def fill(self):
        self.collection = {}
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
                match = expression['sam record'].search(line)
                if match:
                    sam_record = match.groupdict()

                    # Parsing
                    for field in [
                        'FLAG',
                        'POS',
                        'MAPQ',
                        'PNEXT',
                        'TLEN'
                    ]:
                        sam_record[field] = int(sam_record[field])

                    if 'OPT' in sam_record:
                        optional = sam_record['OPT'].strip('\t').split('\t')
                        for o in optional:
                            m = expression['sam optional field'].search(o)
                            if m:
                                o = m.groupdict()
                                if expression['sam optional field value'][o['TYPE']].match(o['VALUE']):
                                    if o['TYPE'] is 'i':
                                        o['VALUE'] = int(o['VALUE'])
                                    elif o['TYPE'] is 'f':
                                        o['VALUE'] = float(o['VALUE'])
                                    sam_record[o['TAG']] = o['VALUE']
                                else:
                                    self.log.error('ignoring invalid %s optional field %s', o['TYPE'], o['VALUE'])

                    if sam_record['QNAME'] in self.collection:
                        record = self.collection[sam_record['QNAME']]
                    else:
                        record = {
                            'OG': 'UK',                         # Organism
                            'CC': 1,                            # Correct Classification
                            'DX': 0,                            # Demultiplexer Exchange
                            'WL': 0,                            # Wrong Lane
                            'WO': 0,                            # Wrong Organism
                            'AS': 0,                            # Organism Alignment Score
                            'MQ': 0,                            # MAPQ
                            'EQ': 0,                            # 0 both agree | 1 don't agree
                            'OBS': sam_record['BC'],            # Observed Barcode Sequence
                            'OBQ': sam_record['QT'],            # Observed Barcode Quality
                            'BEE': 0,                           # Barcode Expected Error
                            'QS': 0,                            # probabilistic class
                            'PBS': sam_record['XL'],            # Probabilistic Barcode Sequence
                            'BS': sam_record['XL'],             # Binned Barcode Sequence
                            'PBD': 0,                           # Probabilistic Barcode Distance
                            'BD': 0,                            # Binned Barcode Distance
                            'DBS': sam_record['XM'],            # Deterministic Barcode Sequence
                            'DBD': 0,                           # Deterministic Barcode Distance
                            'P': 0,                             # Probability
                            'CP': 0,                            # Conditioned Probability
                            'REE': sam_record['EE'],            # Read Expected Error
                            'RS': sam_record['SEQ'],            # Read Sequence
                            'RQ': sam_record['QUAL'],           # Read Quality
                            'ID': sam_record['QNAME'],          # Read ID
                            'LN': None,                         # Lane
                            'TX': None,                         # Taxonomy ID
                            'TD': None,                         # Taxonomy division ID
                            'SAM': [],
                        }

                        if 'XP' in sam_record:
                            record['CP'] = sam_record['XP']

                        if 'DQ' in sam_record:
                            record['P'] = sam_record['DQ']

                        if 'XC' in sam_record:
                            record['QS'] = sam_record['XC']

                        if 'XD' in sam_record:
                            record['PBD'] = sam_record['XD']
                            record['DBD'] = sam_record['XD']

                        if 'YD' in sam_record:
                            record['BD'] = sam_record['YD']

                        record['LN'] = int(record['ID'].split(':')[3])

                        if record['DBS'] != record['PBS']:
                            record['EQ'] = 1

                        barcode = self.barcode[record['PBS']]
                        if any(barcode['LN']) and record['LN'] not in barcode['LN']:
                            record['WL'] = 1

                        # The probabilistic and deterministic do not agree
                        if record['EQ']: 

                            # REMOVED: probabilistic moved the read from a library to undetermined
                            if '=' in record['PBS']:
                                record['DX'] = 1
                                record['BS'] = record['DBS']
                                record['BD'] = record['DBD']

                            # RECOVERED: probabilistic moved the read from a undetermined to a library
                            elif '=' in record['DBS']:
                                record['DX'] = 2

                            # RECLASSIFIED: probabilistic moved the read from a one library to another
                            else:
                                record['DX'] = 3

                        self.collection[record['ID']] = record

                        if self.full:
                            break
                    record['SAM'].append(line)
                else:
                    self.log.error('invalid sam syntax %s', line)

        return not self.empty

    def to_dusted_fasta(self, resolution):
        buffer = StringIO();
        position = 0
        read_id = None
        seq = None
        for record in self.collection.values():
            dusted = record['RS'] if 'DS' not in record else record['DS']
            if not position:
                read_id = record['ID']
                seq = dusted
            else:
                seq += 'N'
                seq += dusted

            position = (position + 1) % resolution
            if not position:
                buffer.write('>{}\n{}\n'.format(read_id, seq))
                read_id = None
                seq = None
        buffer.seek(0)
        return buffer

    def to_fastq(self):
        buffer = StringIO();
        for record in self.collection.values():
            buffer.write('@{}\n{}\n+\n{}\n'.format(record['ID'], record['RS'], record['RQ']))
        buffer.seek(0)
        return buffer

    def to_dusted_fastq(self):
        buffer = StringIO();
        for record in self.collection.values():
            if record['TX'] is None:
                # only align records that kraken hasn't already assigned
                dusted = record['RS'] if 'DS' not in record else record['DS']
                buffer.write('@{}\n{}\n+\n{}\n'.format(record['ID'], dusted, record['RQ']))
        buffer.seek(0)
        return buffer

    def to_sam(self):
        buffer = StringIO();
        for record in self.collection.values():
            for sam_record in record['SAM']:
                buffer.write(sam_record)
        buffer.seek(0)
        return buffer

    def to_fasta(self):
        buffer = StringIO();
        for record in self.collection.values():
            buffer.write('>{}\n{}\n'.format(record['ID'], record['RS']))
        buffer.seek(0)
        return buffer

    def bwa(self, name):
        threads = int(1.25 * THREADS)
        o = self.configuration['organism'][name]
        command = [
            BWA, 'mem',
            '-k', '10',             # minimum seed length
            '-T', '24',             # minimum score to output
            '-t', str(threads),     # number of threads
            '-B', '2',              # penalty for a mismatch
            o['path'],              # organism reference fasta
            '-',                    # output
        ]

        process = Popen(
            args=command,
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )

        # fastq
        buffer = self.to_dusted_fastq()
        output, error = process.communicate(input=buffer.read().encode('utf8'))
        if output:
            buffer = StringIO(output.decode('utf8'))
            for line in buffer:
                hit = None
                if line and line[0] != '@':
                    line = line.split()
                    # 1     QNAME   String
                    # 2     FLAG    Int
                    # 3     RNAME   String
                    # 4     POS     Int
                    # 5     MAPQ    Int
                    # 6     CIGAR   String
                    # 7     RNEXT   String
                    # 8     PNEXT   Int
                    # 9     TLEN    Int
                    # 10    SEQ     String
                    # 11    QUAL    String
                    try:
                        flag = int(line[1])
                        if not flag & 0x4:
                            hit = {
                                'ID': line[0],
                                'FLAG': flag,
                                'MAPQ': int(line[4]),
                            }
                            record = self.collection[hit['ID']]
                            optional = dict([(v[0],v) for v in [i.split(':') for i in line[11:]]])
                            if 'AS' in optional:
                                hit['AS'] = int(optional['AS'][2])

                            if 'NM' in optional:
                                hit['NM'] = int(optional['NM'][2])
                    except ValueError:
                        pass

                    if hit and 'AS' in hit and hit['AS'] > record['AS']:
                        # if record['OG'] != 'UK':
                        #     print('SWITCH {} with {} to {} with {}'.format(record['OG'], record['AS'], name, hit['AS']))

                        record['OG'] = name
                        record['AS'] = hit['AS']
                        record['MQ'] = hit['MAPQ']

    def dust(self):
        def update(d):
            if d['DS'] and self.lowComplexity.search(d['DS']):
                record = self.collection[d['ID']]
                record['DS'] = self.lowComplexity.sub('N', d['DS'])

        command = [
            DUSTMASKER,
            '-outfmt', 'fasta',
            '-level', '20'
        ]

        process = Popen(
            args=command,
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )

        output, error = process.communicate(input=self.to_fasta().read().encode('utf8'))
        if output:
            buffer = StringIO(output.decode('utf8'))
            d = None
            for line in buffer:
                line = line.strip()
                if line:
                    if line[0] == '>':
                        if d is not None:
                            update(d)
                        d = { 'ID': line[1:], 'DS': '' }

                    elif d is not None:
                        d['DS'] += line

            if d is not None:
                update(d)

    def kraken(self):
        # /dev/fd/0 is standard input or file descriptor 0
        command = [
            KRAKEN,
            '--db', KRAKEN_DB,
            '--only-classified-output',
            '--threads', str(THREADS),
            '/dev/fd/0'
        ]
        process = Popen(
            args=command,
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate(input=self.to_dusted_fasta(2).read().encode('utf8'))
        if output:
            buffer = StringIO(output.decode('utf8'))
            record = None
            for line in buffer:
                line = line.strip('\n').strip()
                if line:
                    k_record = line.split()
                    read_id = k_record[1]
                    taxonomy_id = int(k_record[2])
                    if taxonomy_id in self.taxonomy:
                        division_id = self.taxonomy[taxonomy_id]
                        record = self.collection[read_id]
                        if division_id != 8 and division_id in self.division:
                            # if this is not unassigned

                            record['TX'] = taxonomy_id
                            record['TD'] = division_id

                            if taxonomy_id == 374840:
                                # is this PhiX
                                record['OG'] = 'PX'
                            else:
                                record['OG'] = self.division[division_id]['key']
                            # print('\t'.join([str(i) for i in [read_id, record['TX'], record['TD'], record['OG'] ]]))

    def filtered_kraken(self):
        command = [
            KRAKEN,
            '--db', KRAKEN_DB,
            '--only-classified-output',
            '--threads', str(THREADS),
            '/dev/fd/0'
        ]
        process = Popen(
            args=command,
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate(input=self.to_dusted_fasta(2).read().encode('utf8'))

        command = [
            KRAKEN_FILTER,
            '--db', KRAKEN_DB,
            # '--threshold', '0.4'
        ]
        process = Popen(
            args=command,
            env=None,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate(input=output)
        if output:
            buffer = StringIO(output.decode('utf8'))
            record = None
            for line in buffer:
                line = line.strip('\n').strip()
                if line:
                    k_record = line.split()
                    read_id = k_record[1]
                    taxonomy_id = int(k_record[2])
                    if taxonomy_id in self.taxonomy:
                        division_id = self.taxonomy[taxonomy_id]
                        record = self.collection[read_id]
                        if division_id != 8:
                            # if this is not unassigned

                            record['TX'] = taxonomy_id
                            record['TD'] = division_id
                            record['KP'] = float(k_record[4].strip('P='))   # Kraken P

                            if taxonomy_id == 374840:
                                # is this PhiX
                                record['OG'] = 'PX'

                            elif division_id in self.division:
                                record['OG'] = self.division[division_id]['key']
                            # print('\t'.join([str(i) for i in [read_id, record['TX'], record['TD'], record['KP'], record['OG'] ]]))


def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    flowcell = configuration['flowcell'][FLOWCELL]

    if(len(sys.argv)>1):
        name = sys.argv[1]
        if name in configuration['flowcell']:
            flowcell = configuration['flowcell'][name]
            if(len(sys.argv)>2):
                buffer_size = sys.argv[2]
                try:
                    flowcell['capacity'] = int(buffer_size)
                except ValueError as e:
                    print('{} is not a valid buffer size'.format(buffer_size))
                    exit()
        else:
            print('unknown flowcell {}'.format(name))
            exit()
    pipeline = Pipeline(configuration, flowcell)
    pipeline.run()

if __name__ == '__main__':
    main()
