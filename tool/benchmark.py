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

import io
import os
import re
import sys
import json
import numpy
import pickle
import logging
import hashlib
from copy import deepcopy
from datetime import datetime, date

from core.error import *
from core import log_levels
from core import CommandLineParser
from core import Job
from core import merge
from core import to_json
from core import prepare_path

from simulation import EstimatePrior, AdjustPrior
from simulation import SimulateBarcode
from simulation import SimulateSubstitution
from simulation import PamldDemultiplex
from simulation import PamldAccuratePriorDemultiplex
from simulation import PamldUniformDemultiplex
from simulation import MddDemultiplex
from simulation import DemlDemultiplex
from simulation import ToDeML
from simulation import Analyze
from simulation import Summarize

class Benchmark(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.instruction['database path'] = os.path.join(self.home, 'benchmark.json')
        self.instruction['pickeled database path'] = os.path.join(self.home, 'benchmark.pickle')
        self.ontology['db'] = None
        self.ontology['db sha1'] = None

    @property
    def db(self):
        return self.ontology['db']

    def load_db(self):
        if os.path.exists(self.instruction['pickeled database path']):
            self.log.debug('loading pickeled database %s', self.instruction['pickeled database path'])
            with io.open(self.instruction['pickeled database path'], 'rb') as file:
                pickled = pickle.load(file)
                self.ontology['db'] = pickled['db']
                self.ontology['db sha1'] = pickled['db sha1']

        if os.path.exists(self.instruction['database path']):
            self.log.debug('reading json database %s', self.instruction['database path'])
            with io.open(self.instruction['database path'], 'rb') as file:
                content = file.read()
                checksum = hashlib.sha1(content).hexdigest()
                if self.ontology['db sha1'] != checksum:
                    try:
                        self.log.debug('decoding json database %s', self.instruction['database path'])
                        self.ontology['db'] = json.loads(content.decode('utf8'))
                        self.ontology['db sha1'] = checksum
                    except json.decoder.JSONDecodeError as e:
                        self.log.warning('ignoring corrupt database %s', self.instruction['database path'])
                else:
                    self.log.debug('skipping json database decoding %s', self.instruction['database path'])

        if self.db is None:
            self.ontology['db'] = { 'created': str(datetime.now()) }
            self.ontology['db sha1'] = None

        if 'simulation' not in self.db:
            self.db['simulation'] = {}

    def save_db(self, clear=False):
        def persist_pickle(checksum):
            prepare_path(self.instruction['pickeled database path'], self.log)
            with io.open(self.instruction['pickeled database path'], 'wb') as file:
                self.log.info('persisting pickeled database %s', self.instruction['pickeled database path'])
                pickled = {
                    'db': self.db,
                    'db sha1': checksum,
                }
                pickle.dump(pickled, file)

        content = to_json(self.db).encode('utf8')
        checksum  = hashlib.sha1(content).hexdigest()
        if checksum != self.ontology['db sha1']:
            prepare_path(self.instruction['database path'], self.log)
            with io.open(self.instruction['database path'], 'wb') as file:
                self.log.info('persisting database %s', self.instruction['database path'])
                file.write(content)
            self.ontology['db sha1'] = checksum
            persist_pickle(checksum)
        else:
            self.log.debug('skipping json db flush')

        if not os.path.exists(self.instruction['pickeled database path']):
            persist_pickle(checksum)

    def execute(self):
        self.load_db()
        if self.action == 'plan':
            self.plan(self.ontology)

        if self.action == 'rebuild':
            self.rebuild_db(self.ontology)

        if self.action == 'simulate_barcode':
            self.simulate_barcode(self.ontology)

        elif self.action == 'simulate_substitution':
            self.simulate_substitution(self.ontology)

        elif self.action == 'todeml':
            self.todeml(self.ontology)

        elif self.action == 'sense_prior':
            self.sense_prior(self.ontology)

        elif self.action == 'adjust_prior':
            self.adjust_prior(self.ontology)

        elif self.action == 'demux_pamld':
            self.demux_pamld(self.ontology)

        elif self.action == 'demux_mdd':
            self.demux_mdd(self.ontology)

        elif self.action == 'demux_deml':
            self.demux_deml(self.ontology)

        elif self.action == 'analyze_pamld':
            self.analyze_pamld(self.ontology)

        elif self.action == 'analyze_mdd':
            self.analyze_mdd(self.ontology)

        elif self.action == 'analyze_deml':
            self.analyze_deml(self.ontology)

        elif self.action == 'collect':
            self.collect(self.ontology)

        elif self.action == 'summarize':
            self.summarize(self.ontology)
        # self.save_db()

    def plan(self, ontology):
        path = os.path.realpath(self.instruction['path'])
        if os.path.exists(path):
            self.log.debug('loading simulation plan %s', self.instruction['path'])
            with io.open(path, 'rb') as file:
                plan = json.loads(file.read().decode('utf8'))

                if 'job' in plan:
                    compiled = []
                    for job in plan['job']:
                        instruction = merge(self.instruction, job)
                        if 'default' in plan:
                            instruction = merge(plan['default'], instruction)
                        compiled.append({ 'instruction': instruction })

                    for job in compiled:
                        if 'action' in job['instruction']:
                            if job['instruction']['action'] == 'benchmark_substitution':
                                self.benchmark_substitution(job)

                            elif job['instruction']['action'] == 'simulate_barcode':
                                self.simulate_barcode(job)

                            elif job['instruction']['action'] == 'simulate_substitution':
                                self.simulate_substitution(job)

                            elif job['instruction']['action'] == 'todeml':
                                self.todeml(self.ontology)

                            elif job['instruction']['action'] == 'sense_prior':
                                self.sense_prior(job)

                            elif job['instruction']['action'] == 'adjust_prior':
                                self.adjust_prior(job)

                            elif job['instruction']['action'] == 'demux_pamld':
                                self.demux_pamld(job)

                            elif job['instruction']['action'] == 'demux_acurate_prior_pamld':
                                self.demux_acurate_prior_pamld(job)

                            elif job['instruction']['action'] == 'demux_uniform_prior_pamld':
                                self.demux_uniform_prior_pamld(job)

                            elif job['instruction']['action'] == 'demux_mdd':
                                self.demux_mdd(job)

                            elif job['instruction']['action'] == 'demux_deml':
                                self.demux_deml(job)

                            elif job['instruction']['action'] == 'analyze_pamld':
                                self.analyze_pamld(job)

                            elif job['instruction']['action'] == 'analyze_mdd':
                                self.analyze_mdd(job)

                            elif job['instruction']['action'] == 'analyze_deml':
                                self.analyze_deml(job)

        else: raise NoConfigurationFileError('plan file {} not found'.format(path))

    def benchmark_substitution(self, ontology):
        o = deepcopy(ontology)
        o['instruction']['action'] = 'simulate_substitution'
        job = self.simulate_substitution(o)
        ssid = job.ssid

        o = deepcopy(ontology)
        o['instruction']['action'] = 'todeml'
        o['instruction']['ssid'] = ssid
        job = self.todeml(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'sense_prior'
        o['instruction']['ssid'] = ssid
        job = self.sense_prior(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'adjust_prior'
        o['instruction']['ssid'] = ssid
        job = self.adjust_prior(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'demux_pamld'
        o['instruction']['ssid'] = ssid
        job = self.demux_pamld(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'analyze_pamld'
        o['instruction']['ssid'] = ssid
        job = self.analyze_pamld(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'demux_acurate_prior_pamld'
        o['instruction']['ssid'] = ssid
        job = self.demux_acurate_prior_pamld(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'analyze_acurate_prior_pamld'
        o['instruction']['ssid'] = ssid
        job = self.analyze_acurate_prior_pamld(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'demux_uniform_prior_pamld'
        o['instruction']['ssid'] = ssid
        job = self.demux_uniform_prior_pamld(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'analyze_uniform_pamld'
        o['instruction']['ssid'] = ssid
        job = self.analyze_uniform_pamld(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'demux_mdd'
        o['instruction']['ssid'] = ssid
        job = self.demux_mdd(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'analyze_mdd'
        o['instruction']['ssid'] = ssid
        job = self.analyze_mdd(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'demux_deml'
        o['instruction']['ssid'] = ssid
        job = self.demux_deml(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'analyze_deml'
        o['instruction']['ssid'] = ssid
        job = self.analyze_deml(o)

    def simulate_barcode(self, ontology):
        job = None
        self.log.info('simulating barcode indices')

        if 'bsid' in ontology['instruction'] and ontology['instruction']['bsid'] in self.db['simulation']:
            barcode_simulation_node = self.db['simulation'][ontology['instruction']['bsid']]
            ontology['model'] = deepcopy(barcode_simulation_node['barcode']['model'])

        job = SimulateBarcode(ontology)
        job.execute()

        if job.is_model_dirty:
            if job.bsid not in self.db['simulation']:
                self.db['simulation'][job.bsid] = {}

            if 'barcode' not in self.db['simulation'][job.bsid]:
                self.db['simulation'][job.bsid]['barcode'] = {}

            if 'substitution' not in self.db['simulation'][job.bsid]:
                self.db['simulation'][job.bsid]['substitution'] = {}

            self.db['simulation'][job.bsid]['barcode']['model'] = job.model_summary

            self.save_db()
        return job

    def simulate_substitution(self, ontology):
        job = None
        self.log.info('simulating errors on barcode indices')

        # load the barcode simulation model
        if 'bsid' in ontology['instruction'] and ontology['instruction']['bsid'] in self.db['simulation']:
            barcode_simulation_node = self.db['simulation'][ontology['instruction']['bsid']]
            if 'ssid' in ontology['instruction'] and ontology['instruction']['ssid'] in barcode_simulation_node['substitution']:
                substitution_simulation_node = barcode_simulation_node['substitution'][ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(substitution_simulation_node['model'])
            else:
                ontology['model'] = deepcopy(barcode_simulation_node['barcode']['model'])

            job = SimulateSubstitution(ontology)
            job.execute()
            if job.is_model_dirty:
                barcode_simulation_node['substitution'][job.ssid] = {
                    'model': job.model_summary,
                }
                self.save_db()
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def todeml(self, ontology):
        job = None
        self.log.info('transcoding simulated data to deML syntax')

        # load the barcode simulation model
        if 'bsid' in ontology['instruction'] and ontology['instruction']['bsid'] in self.db['simulation']:
            barcode_simulation_node = self.db['simulation'][ontology['instruction']['bsid']]
            if 'ssid' in ontology['instruction'] and ontology['instruction']['ssid'] in barcode_simulation_node['substitution']:
                substitution_simulation_node = barcode_simulation_node['substitution'][ontology['instruction']['ssid']]
                if 'model' in substitution_simulation_node:
                    ontology['model'] = deepcopy(substitution_simulation_node['model'])
                    job = ToDeML(ontology)
                    job.execute()
        return job

    def sense_prior(self, ontology):
        job = None
        self.log.info('estimating priors')

        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = EstimatePrior(ontology)
                job.execute()

                is_dirty = False
                if job.is_model_dirty:
                    experiment['model'] = job.model_summary
                    is_dirty = True

                if job.is_execution_summary_dirty:
                    experiment['sense prior execution'] = job.execution_summary
                    is_dirty = True

                if is_dirty:
                    self.save_db()
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def adjust_prior(self, ontology):
        job = None
        self.log.info('adjusting priors')

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = AdjustPrior(ontology)
                job.execute()
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])
        return job

    def demux_pamld(self, ontology):
        job = None
        self.log.info('demultiplexing pamld')

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = PamldDemultiplex(ontology)
                job.execute()
                if job.is_execution_summary_dirty:
                    experiment['pamld demultiplex execution'] = job.execution_summary
                    self.save_db()
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def demux_acurate_prior_pamld(self, ontology):
        job = None
        self.log.info('demultiplexing pamld with accurate prior')

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = PamldAccuratePriorDemultiplex(ontology)
                job.execute()
                if job.is_execution_summary_dirty:
                    experiment['pamld accurate prior demultiplex execution'] = job.execution_summary
                    self.save_db()
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def demux_uniform_prior_pamld(self, ontology):
        job = None
        self.log.info('demultiplexing pamld with uniform prior')

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = PamldUniformDemultiplex(ontology)
                job.execute()
                if job.is_execution_summary_dirty:
                    experiment['pamld uniform demultiplex execution'] = job.execution_summary
                    self.save_db()
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def demux_mdd(self, ontology):
        job = None
        self.log.info('demultiplexing mdd')

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = MddDemultiplex(ontology)
                job.execute()
                if job.is_execution_summary_dirty:
                    experiment['mdd demultiplex execution'] = job.execution_summary
                    self.save_db()
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def demux_deml(self, ontology):
        job = None
        self.log.info('demultiplexing with deml')

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = DemlDemultiplex(ontology)
                job.execute()
                if job.is_execution_summary_dirty:
                    experiment['deml demultiplex execution'] = job.execution_summary
                    self.save_db()
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def analyze_pamld(self, ontology):
        job = None
        self.log.info('analyzing pamld results')
        ontology['instruction']['tool'] = 'pamld'

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            barcode_simulation = self.db['simulation'][ontology['instruction']['bsid']]
            if ontology['instruction']['ssid'] in barcode_simulation['substitution']:
                substitution_simulation = barcode_simulation['substitution'][ontology['instruction']['ssid']]
                if 'pamld demultiplex analysis' not in substitution_simulation:
                    ontology['model'] = deepcopy(substitution_simulation['model'])
                    ontology['instruction']['input'] = ontology['model']['location']['pamld demultiplex path']

                    job = Analyze(ontology)
                    job.execute()
                    substitution_simulation['pamld demultiplex analysis'] = job.model_summary
                    self.save_db()
                else:
                    self.log.info('skipping pamld analysis for %s', ontology['instruction']['ssid'])
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def analyze_acurate_prior_pamld(self, ontology):
        job = None
        self.log.info('analyzing acurate prior pamld results')
        ontology['instruction']['tool'] = 'pamld_ap'

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            barcode_simulation = self.db['simulation'][ontology['instruction']['bsid']]
            if ontology['instruction']['ssid'] in barcode_simulation['substitution']:
                substitution_simulation = barcode_simulation['substitution'][ontology['instruction']['ssid']]
                if 'pamld accurate prior demultiplex analysis' not in substitution_simulation:
                    ontology['model'] = deepcopy(substitution_simulation['model'])
                    ontology['instruction']['input'] = ontology['model']['location']['pamld accurate prior demultiplex path']

                    job = Analyze(ontology)
                    job.execute()
                    substitution_simulation['pamld accurate prior demultiplex analysis'] = job.model_summary
                    self.save_db()
                else:
                    self.log.info('skipping pamld accurate prior analysis for %s', ontology['instruction']['ssid'])
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def analyze_uniform_pamld(self, ontology):
        job = None
        self.log.info('analyzing uniform pamld results')
        ontology['instruction']['tool'] = 'pamld_u'

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            barcode_simulation = self.db['simulation'][ontology['instruction']['bsid']]
            if ontology['instruction']['ssid'] in barcode_simulation['substitution']:
                substitution_simulation = barcode_simulation['substitution'][ontology['instruction']['ssid']]
                if 'pamld uniform demultiplex analysis' not in substitution_simulation:
                    ontology['model'] = deepcopy(substitution_simulation['model'])
                    ontology['instruction']['input'] = ontology['model']['location']['pamld uniform demultiplex path']

                    job = Analyze(ontology)
                    job.execute()
                    substitution_simulation['pamld uniform demultiplex analysis'] = job.model_summary
                    self.save_db()
                else:
                    self.log.info('skipping pamld uniform analysis for %s', ontology['instruction']['ssid'])
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def analyze_mdd(self, ontology):
        job = None
        self.log.info('analyzing mdd results')
        ontology['instruction']['tool'] = 'mdd'

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            barcode_simulation = self.db['simulation'][ontology['instruction']['bsid']]
            if ontology['instruction']['ssid'] in barcode_simulation['substitution']:
                substitution_simulation = barcode_simulation['substitution'][ontology['instruction']['ssid']]
                if 'mdd demultiplex analysis' not in substitution_simulation:
                    ontology['model'] = deepcopy(substitution_simulation['model'])
                    ontology['instruction']['input'] = ontology['model']['location']['mdd demultiplex path']
                    job = Analyze(ontology)
                    job.execute()
                    substitution_simulation['mdd demultiplex analysis'] = job.model_summary
                    self.save_db()
                else:
                    self.log.info('skipping mdd analysis for %s', ontology['instruction']['ssid'])
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def analyze_deml(self, ontology):
        job = None
        self.log.info('analyzing deML results')
        ontology['instruction']['tool'] = 'deml'

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            barcode_simulation = self.db['simulation'][ontology['instruction']['bsid']]
            if ontology['instruction']['ssid'] in barcode_simulation['substitution']:
                substitution_simulation = barcode_simulation['substitution'][ontology['instruction']['ssid']]
                if 'deml demultiplex analysis' not in substitution_simulation:
                    ontology['model'] = deepcopy(substitution_simulation['model'])
                    ontology['instruction']['input'] = ontology['model']['location']['deml demultiplex path']
                    job = Analyze(ontology)
                    job.execute()
                    substitution_simulation['deml demultiplex analysis'] = job.model_summary
                    self.save_db()
                else:
                    self.log.info('skipping deml analysis for %s', ontology['instruction']['bsid'])
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        return job

    def summarize(self, ontology):
        job = None
        if self.instruction['bsid'] in self.db['simulation']:
            ontology['barcode simulation'] = self.db['simulation'][self.instruction['bsid']]
            job = Summarize(ontology)
            job.execute()
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        self.save_db()
        return job

    def rebuild_db(self, ontology):
        def compile_model(node):
            default = {
                'genealogy': {},
                'location': {},
            }
            model = merge(default, node)

            if 'bsid' in model['genealogy']:
                model['location']['barcode simulation home'] = os.path.join('simulation', model['genealogy']['bsid'])
                for key, relative in {
                    'deml index path': 'deml_index.txt',
                    'mdd configuration path': 'mdd_configuration.json',
                    'pamld accurate prior configuration path': 'pamld_accurate_prior_configuration.json',
                    'pamld uniform configuration path': 'pamld_uniform_configuration.json',
                    'simulated barcode path': 'simulated_barcode.bam',
                }.items():
                    if key not in model['location']:
                        model['location'][key] = os.path.join(model['location']['barcode simulation home'], relative)

                if 'ssid' in model['genealogy']:
                    model['location']['substitution simulation home'] = os.path.join(model['location']['barcode simulation home'], 'substitution', model['genealogy']['ssid'])
                    for key, relative in {
                        'deml demultiplex path' : 'deml_demultiplex.bam',
                        'deml simulated substitution path' : 'deml_simulated_substitution.bam',
                        'deml summary path' : 'deml_summary.txt',
                        'mdd demultiplex path' : 'mdd_demultiplex.bam',
                        'mdd demultiplex report path' : 'mdd_demultiplex_report.json',
                        'pamld adjusted configuration path' : 'pamld_adjusted_configuration.json',
                        'pamld demultiplex path' : 'pamld_demultiplex.bam',
                        'pamld uniform demultiplex path' : 'pamld_uniform_demultiplex.bam',
                        'pamld uniform demultiplex report path' : 'pamld_uniform_demultiplex_report.json',
                        'pamld demultiplex report path' : 'pamld_demultiplex_report.json',
                        'pamld prior estimate path' : 'prior_estimate.json',
                        'pamld accurate prior demultiplex path' : 'pamld_accurate_prior_demultiplex.bam',
                        'pamld accurate prior demultiplex report path' : 'pamld_accurate_prior_demultiplex_report.json',
                        'simulated substitution path' : 'simulated_substitution.bam',
                    }.items():
                        if key not in model['location']:
                            model['location'][key] = os.path.join(model['location']['substitution simulation home'], relative)

                    # for key in [
                    # ]:
                    #     if key in model['location']:
                    #         del model['location'][key]

            return model

        def compile_barcode_simulation(bsid):
            if bsid in self.db['simulation']:
                barcode_simulation_node = self.db['simulation'][bsid]
                barcode_simulation_node['barcode']['model'] = compile_model(barcode_simulation_node['barcode']['model'])

                if 'substitution' in barcode_simulation_node:
                    self.log.debug('compiling %d substitution in simulation %s', len(barcode_simulation_node['substitution']), bsid)
                    for ssid, substitution_simulation_node in barcode_simulation_node['substitution'].items():
                        substitution_simulation_node['model'] = compile_model(substitution_simulation_node['model'])
                        self.log.debug('compiling substitution %s', ssid)

                        for key in [
                            'deml demultiplex analysis',
                            'mdd demultiplex analysis',
                            'pamld demultiplex analysis',
                            'pamld uniform demultiplex analysis',
                            'pamld accurate prior demultiplex analysis',
                        ]:
                            if key in substitution_simulation_node:
                                self.log.debug('compiling %s', key)
                                substitution_simulation_node[key] = compile_model(substitution_simulation_node[key])

        self.log.info('cleaning simulation database')
        if 'bsid' in ontology['instruction']:
            compile_barcode_simulation(ontology['instruction']['bsid'])
        else:
            for bsid in self.db['simulation'].keys():
                compile_barcode_simulation(bsid)

        self.save_db()

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    job = None

    try:
        command = CommandLineParser('benchmark')
        if command.help_triggered:
            command.help()
            sys.exit(0)
        else:
            if 'verbosity' in command.instruction and command.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

            job = Benchmark(command.configuration)
            job.execute()

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
        if job: job.close()

    sys.exit(0)

if __name__ == '__main__':
    main()
