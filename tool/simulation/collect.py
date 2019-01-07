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

class Collect(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.analysis_path_re = re.compile(self.instruction['pattern'])
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
            self.log.debug('creating a new database %s', self.db_path)
            self.ontology['db'] = { 'experiment': {}, 'created': str(datetime.now()) }

        self.ontology['db']['loaded'] = str(datetime.now())

    def persist(self):
        prepare_path(self.db_path, self.log)
        with io.open(self.db_path, 'wb') as file:
            self.log.debug('persisting database %s', self.db_path)
            self.db['saved'] = str(datetime.now())
            content = json.dumps(self.db, sort_keys=True, ensure_ascii=False, indent=4).encode('utf8')
            file.write(content)

    def update(self, experiment):
        id = experiment['genealogy']['error simulation id']
        experiment['updated'] = str(datetime.now())
        if id not in self.db['experiment']:
            self.db['experiment'][id] = []
        self.db['experiment'][id].append(experiment)

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

    def collect_analysis(self):
        dirname = self.instruction['path']
        for basename in os.listdir(dirname):
            path = os.path.abspath(os.path.join(dirname, basename))
            if inode_type(path) == 'file':
                match = self.analysis_path_re.search(basename)
                if match:
                    experiment = {
                        'basename': basename,
                        'path': path,
                    }
                    self.log.info('collecting results from %s', experiment['basename'])
                    with io.open(experiment['path'], 'rb') as file:
                        node = json.loads(file.read().decode('utf8'))

                        if 'genealogy' in node and 'error simulation id' in node['genealogy']:
                            experiment['genealogy'] = deepcopy(node['genealogy'])
                            for topic in [ 'multiplex', 'cellular', 'molecular' ]:
                                if topic in node:

                                    if isinstance(node[topic], dict):
                                        experiment[topic] = self.collect_decoder_analysis(node[topic])

                                    elif isinstance(node[topic], list):
                                        experiment[topic] = []
                                        for decoder in node[topic]:
                                            experiment[topic].append(self.collect_decoder_analysis(decoder))

                            self.update(experiment)
                    # print(to_json(self.ontology['experiment']))
                    # exit(0)

    def execute(self):
        self.load()
        self.collect_analysis()
        self.finalize()
        self.persist()

    def finalize(self):
        pass
