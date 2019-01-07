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

class Summarize(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        if 'database' in self.instruction and self.instruction['database']:
            self.ontology['db path'] = self.instruction['database']
        else:
            self.ontology['db path'] = os.path.join(self.instruction['home'], 'benchmark_db.json')
        self.ontology['db path'] = os.path.abspath(os.path.expanduser(os.path.expandvars(self.ontology['db path'])))

    @property
    def db_path(self):
        return self.ontology['db path']

    @property
    def db(self):
        return self.ontology['db']

    def load(self):
        if os.path.exists(self.db_path):
            self.log.debug('loading existing database %s', self.db_path)
            with io.open(self.db_path, 'rb') as file:
                self.ontology['db'] = json.loads(file.read().decode('utf8'))
        else:
            self.log.error('no database found in %s', self.db_path)

    def make_summary(self):
        self.ontology['summary'] = { 'total': {}, 'barcode': {} }
        for record in self.ontology['experiment']:
            rate = record['rate']
            tool = record['tool']
            self.log.info('collecting data from %s', record['basename'])
            with io.open(record['path'], 'rb') as file:
                experiment = json.loads(file.read().decode('utf8'))
                if rate not in self.summary['total']: self.summary['total'][rate] = {}
                if tool not in self.summary['total'][rate]: self.summary['total'][rate][tool] = {}
                for qc in  self.scale['qc']:
                    self.summary['total'][rate][tool][qc] = experiment['decoder']['A5KVK'][qc]

                for k,v in experiment['decoder']['A5KVK']['codec'].items():
                    if k not in self.summary['barcode']: self.summary['barcode'][k] = {}
                    if rate not in self.summary['barcode'][k]: self.summary['barcode'][k][rate] = {}
                    if tool not in self.summary['barcode'][k][rate]: self.summary['barcode'][k][rate][tool] = {}
                    for qc in  self.scale['qc']:
                        self.summary['barcode'][k][rate][tool][qc] = v[qc]

        for rate in self.summary['total'].values():
            for tool in rate.values():
                both = {}
                for qc in tool.values():
                    for k,v in qc.items():
                        if k not in both: both[k] = v
                        else: both[k] += v
                tool['both'] = both

        for barcode in self.summary['barcode'].values():
            for rate in barcode.values():
                for tool in rate.values():
                    both = {}
                    for qc in tool.values():
                        for k,v in qc.items():
                            if k not in both: both[k] = v
                            else: both[k] += v
                    tool['both'] = both

    def write_summary_csv_for_node(self, name, node):
        output_path = '{}_summary.csv'.format(name)
        content = []
        head = [ 'rate' ]
        for qc in self.header['qc'] + ['B']:
            for field in self.header['field']:
                for tool in self.header[ 'tool' ]:
                    head.append('{}:{}:{}'.format(tool, qc, field))
        content.append(','.join(head))

        for rate in sorted(node.keys()):
            row = [ '{}'.format(rate) ]
            for qc in self.scale['qc'] + ['both']:
                for field in self.scale['field']:
                    for tool in self.scale['tool']:
                        row.append(node[rate][tool][qc][field])
            if row:
                content.append(','.join([str(i) for i in row]))

        with io.open(output_path, 'w') as file:
            file.write('\n'.join(content))

    def write_csv_for_each_barcode(self):
        self.write_summary_csv_for_node('total', self.summary['total'])
        for barcode in sorted(self.summary['barcode'].keys()):
            self.write_summary_csv_for_node(barcode,self.summary['barcode'][barcode])

    def write_summary_csv(self):
        output_path = 'summary.csv'
        content = []
        head = [ 'barcode', 'rate' ]
        for qc in self.header['qc']:
            for field in self.header['field']:
                for tool in self.header[ 'tool' ]:
                    head.append('{}:{}:{}'.format(tool, qc, field))
        content.append(','.join(head))

        for barcode in sorted(self.summary['barcode'].keys()):
            for rate in sorted(self.summary['barcode'][barcode].keys()):
                row = [ barcode, '{}'.format(rate) ]
                for qc in self.scale['qc']:
                    for field in self.scale['field']:
                        for tool in self.scale['tool']:
                            row.append(self.summary['barcode'][barcode][rate][tool][qc][field])
                content.append(','.join([str(i) for i in row]))

        with io.open(output_path, 'w') as file:
            file.write('\n'.join(content))

    def execute(self):
        self.load()

        self.write_csv_for_each_barcode()
        self.write_summary_csv()

        self.finalize()

    def collect_per_barcode():
        pass

    def collect_decoder_analysis(self, decoder):
        meta_decoder = {}
        for item in [
            'key',
            'index',
            'noise',
            'macro',
            'qc fail',
            'qc pass',
            'classified count',
            'simulated classified nucleotide',
            'simulated classified substitution',
            'simulated noise',
            'simulated noise deviation',
            'simulated nucleotide',
            'simulated substitution',
            'substitution rate per classified nucleotide',
            'substitution rate per nucleotide',
        ]:
            if item in decoder:
                meta_decoder[item] = decoder[item]

        if 'codec' in decoder:
            meta_decoder['codec'] = {}
            for barcode in decoder['codec'].values():
                meta_barcode = {}
                meta_decoder['codec'][barcode['key']] = meta_barcode
                for item in [
                    'key',
                    'RG'
                    'index',
                    'barcode',
                    'concentration',
                    'count',
                    'qc fail',
                    'qc pass',
                    'simulated concentration',
                    'simulated nucleotide',
                    'simulated substitution',
                    'substitution rate per barcode',
                    'substitution rate per nucleotide'
                ]:
                    if item in barcode:
                        meta_barcode[item] = barcode[item]

        if 'unclassified' in decoder:
            meta_decoder['unclassified'] = {}
            for item in [
                'RG'
                'index',
                'concentration',
                'count',
                'qc fail',
                'qc pass',
                'simulated concentration',
                'simulated nucleotide',
                'simulated substitution',
                'substitution rate per barcode',
                'substitution rate per nucleotide'
            ]:
                if item in decoder['unclassified']:
                    meta_decoder['unclassified'][item] = decoder['unclassified'][item]
        return meta_decoder

    def execute(self):
        self.load()
        self.collect_analysis()
        self.finalize()
        self.persist()

    def finalize(self):
        pass
