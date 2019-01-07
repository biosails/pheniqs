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
import logging
from copy import deepcopy
from datetime import datetime, date

from core.error import *
from core import CommandLineParser
from core import Job
from core import merge
from core import to_json
from core import prepare_path

from simulation import SensePrior, AdjustPrior
from simulation import SimulateBarcode
from simulation import SimulateSubstitution
from simulation import PheniqsDemultiplex
from simulation import DemlDemultiplex
from simulation import ToDeML
from simulation import Analyze
# from simulation import Collect
# from simulation import Summarize

class Benchmark(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.instruction['database path'] = os.path.join(self.home, 'benchmark.json')

    @property
    def db(self):
        return self.ontology['db']

    def load_db(self):
        if os.path.exists(self.instruction['database path']):
            self.log.debug('loading existing database %s', self.instruction['database path'])
            with io.open(self.instruction['database path'], 'rb') as file:
                try:
                    self.ontology['db'] = json.loads(file.read().decode('utf8'))
                except json.decoder.JSONDecodeError as e:
                    self.log.warning('ignoring corrupt database %s', self.instruction['database path'])

        if 'db' not in self.ontology:
            self.ontology['db'] = { 'created': str(datetime.now()) }

        self.db['loaded'] = str(datetime.now())
        if 'simulation' not in self.db:
            self.db['simulation'] = {}

    def persist_db(self):
        prepare_path(self.instruction['database path'], self.log)
        with io.open(self.instruction['database path'], 'wb') as file:
            self.log.debug('persisting database %s', self.instruction['database path'])
            self.db['saved'] = str(datetime.now())
            file.write(to_json(self.db).encode('utf8'))

    def execute(self):
        self.load_db()
        if self.action == 'plan':
            self.plan(self.ontology)

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

        elif self.action == 'demux_pheniqs':
            self.demux_pheniqs(self.ontology)

        elif self.action == 'demux_deml':
            self.demux_deml(self.ontology)

        elif self.action == 'analyze_pheniqs':
            self.analyze_pheniqs(self.ontology)

        elif self.action == 'analyze_deml':
            self.analyze_deml(self.ontology)

        elif self.action == 'collect':
            self.collect(self.ontology)

        elif self.action == 'summarize':
            self.summarize(self.ontology)
        self.persist_db()

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

                            elif job['instruction']['action'] == 'demux_pheniqs':
                                self.demux_pheniqs(job)

                            elif job['instruction']['action'] == 'demux_deml':
                                self.demux_deml(job)

                            elif job['instruction']['action'] == 'analyze_pheniqs':
                                self.analyze_pheniqs(job)

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
        o['instruction']['action'] = 'demux_pheniqs'
        o['instruction']['ssid'] = ssid
        job = self.demux_pheniqs(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'demux_deml'
        o['instruction']['ssid'] = ssid
        job = self.demux_deml(o)

        o = deepcopy(ontology)
        o['instruction']['action'] = 'analyze_pheniqs'
        o['instruction']['ssid'] = ssid
        job = self.analyze_pheniqs(o)

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

        if job.bsid not in self.db['simulation']:
            self.db['simulation'][job.bsid] = {}

        if 'barcode' not in self.db['simulation'][job.bsid]:
            self.db['simulation'][job.bsid]['barcode'] = {}

        if 'substitution' not in self.db['simulation'][job.bsid]:
            self.db['simulation'][job.bsid]['substitution'] = {}

        self.db['simulation'][job.bsid]['barcode']['model'] = job.summary

        self.persist_db()
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

            barcode_simulation_node['substitution'][job.ssid] = {
                'model': job.summary,
            }
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        self.persist_db()
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

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = SensePrior(ontology)
                job.execute()
                experiment['sense prior execution'] = job.execution
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        self.persist_db()
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

    def demux_pheniqs(self, ontology):
        job = None
        self.log.info('demultiplexing with pheniqs')

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = PheniqsDemultiplex(ontology)
                job.execute()
                experiment['pheniqs demultiplex execution'] = job.execution
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        self.persist_db()
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
                experiment['deml demultiplex execution'] = job.execution
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        self.persist_db()
        return job

    def analyze_pheniqs(self, ontology):
        job = None
        self.log.info('analyzing pheniqs results')
        ontology['instruction']['tool'] = 'pheniqs'

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = Analyze(ontology)
                job.execute()
                experiment['pheniqs demultiplex analysis'] = job.summary
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        self.persist_db()
        return job

    def analyze_deml(self, ontology):
        job = None
        self.log.info('analyzing deML results')
        ontology['instruction']['tool'] = 'deml'

        # load the substitution simulation model
        if ontology['instruction']['bsid'] in self.db['simulation']:
            substitution = self.db['simulation'][ontology['instruction']['bsid']]['substitution']
            if ontology['instruction']['ssid'] in substitution:
                experiment = substitution[ontology['instruction']['ssid']]
                ontology['model'] = deepcopy(experiment['model'])
                job = Analyze(ontology)
                job.execute()
                experiment['deml demultiplex analysis'] = job.summary
            else:
                self.log.error('unknown ssid %s', ontology['instruction']['bsid'])
        else:
            self.log.error('unknown bsid %s', ontology['instruction']['bsid'])

        self.persist_db()
        return job

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
