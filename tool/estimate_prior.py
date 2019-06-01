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
import re
import io
import sys
import json
import logging
from copy import deepcopy
from subprocess import Popen, PIPE

from core.error import *
from core import log_levels
from core import CommandLineParser
from core import Job
from core import merge
from core import to_json

class EstimatePrior(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.log = logging.getLogger('EstimatePrior')

    @property
    def original(self):
        return self.ontology['original']

    @property
    def adjusted(self):
        return self.ontology['adjusted']

    def execute(self):
        self.estimate_prior_from_report()

    def estimate_prior_from_report(self):
        def estimate_classifier_prior(m, r):
            if 'average classified confidence' in r:

                # estimate noise prior
                noise_count = 0
                if 'low conditional confidence count' in r:
                    noise_count += r['low conditional confidence count']

                if 'low confidence count' in r:
                    noise_count += (r['low confidence count'] * (1.0 - r['average classified confidence']))

                m['noise'] = noise_count / r['count']
                not_noise = 1.0 - m['noise']

                # estimate each barcode prior
                barcode_report_by_hash = {}
                for barcode in r['classified']:
                    hash = ''.join(barcode['barcode'])
                    barcode_report_by_hash[hash] = barcode

                if 'codec' in m:
                    for bm in m['codec'].values():
                        hash = ''.join(bm['barcode'])
                        if hash in barcode_report_by_hash:
                            br = barcode_report_by_hash[hash]
                            if 'pf pooled classified fraction' in br and br['pf pooled classified fraction'] > 0:
                                bm['concentration'] = not_noise * br['pf pooled classified fraction']
                            else:
                                bm['concentration'] = 0

        self.make_static_configuration()
        self.load_report()
        self.ontology['adjusted configurtion'] = deepcopy(self.ontology['static configurtion'])

        if self.instruction['split_fastq']:
            self.make_compiled_configuration()
            if 'multiplex' in self.ontology['adjusted configurtion']:
                segment_cardinality = self.ontology['compiled configurtion']['multiplex']['segment cardinality']
                if 'undetermined' in self.ontology['adjusted configurtion']['multiplex']:
                    barcode = self.ontology['adjusted configurtion']['multiplex']['undetermined']
                    barcode['output'] = []
                    for segment_index in range(1,segment_cardinality + 1):
                        barcode['output'].append('{}_undetermined_s{:0>2}.fastq.gz'.format(self.instruction['prefix'], segment_index))

                if 'codec' in self.ontology['adjusted configurtion']['multiplex']:
                    for barcode in self.ontology['adjusted configurtion']['multiplex']['codec'].values():
                        barcode['output'] = []
                        hash = ''.join(barcode['barcode'])
                        for segment_index in range(1,segment_cardinality + 1):
                            barcode['output'].append('{}_{}_s{:0>2}.fastq.gz'.format(self.instruction['prefix'], hash, segment_index))

        elif self.instruction['split_bam']:
            if 'multiplex' in self.ontology['adjusted configurtion']:
                if 'undetermined' in self.ontology['adjusted configurtion']['multiplex']:
                    barcode = self.ontology['adjusted configurtion']['multiplex']['undetermined']
                    barcode['output'] = [ '{}_undetermined.bam'.format(self.instruction['prefix'], segment_index) ]

                if 'codec' in self.ontology['adjusted configurtion']['multiplex']:
                    for barcode in self.ontology['adjusted configurtion']['multiplex']['codec'].values():
                        hash = ''.join(barcode['barcode'])
                        barcode['output'] = [ '{}_{}.bam'.format(self.instruction['prefix'], hash) ]

        for classifier_type in [ 'multiplex', 'cellular', 'molecular' ]:
            if classifier_type in self.ontology['estimation report'] and classifier_type in self.ontology['adjusted configurtion']:
                m = self.ontology['adjusted configurtion'][classifier_type]
                r = self.ontology['estimation report'][classifier_type]

                if isinstance(m, dict):
                    estimate_classifier_prior(m, r)

                elif isinstance(m, list):
                    index = 0
                    model_by_index = {}
                    for item in m:
                        item['index'] = index
                        model_by_index[index] = item
                        index += 1

                    for report_item in r:
                        model_item = model_by_index[report_item['index']]
                        estimate_classifier_prior(model_item, report_item)

        print(to_json(self.ontology['adjusted configurtion']))

    def make_static_configuration(self):
        self.location['original configuration'] = os.path.realpath(os.path.abspath(os.path.expanduser(os.path.expandvars(self.instruction['configuration']))))
        if os.path.exists(self.location['original configuration']):
            command = [ 'pheniqs', 'mux', '--static' ]
            command.append('--config')
            command.append(self.location['original configuration'])

            if self.instruction['sense_input']:
                command.append('--sense-input')

            if 'base_input' in self.instruction:
                command.append('--base-input')
                command.append(self.instruction['base_input'])

            if 'base_output' in self.instruction:
                command.append('--base-output')
                command.append(self.instruction['base_output'])

            if 'input' in self.instruction:
                for value in self.instruction['input']:
                    command.append('--input')
                    command.append(value)

            process = Popen(
                args=command,
                cwd=self.current_working_directoy,
                stdout=PIPE,
                stderr=PIPE
            )
            output, error = process.communicate()
            return_code = process.returncode

            if return_code == 0:
                compiled = json.loads(output.decode('utf8'))
                self.ontology['static configurtion'] = compiled
            else:
                raise BadConfigurationError('pheniqs returned {} when creating a static configuration'.format(return_code))
                for line in error.decode('utf8').splitlines():
                    print(error.decode('utf8'))
        else:
            raise BadConfigurationError('could not find configuration file {}'.format(self.instruction['configurtion']))

    def make_compiled_configuration(self):
        self.location['original configuration'] = os.path.realpath(os.path.abspath(os.path.expanduser(os.path.expandvars(self.instruction['configuration']))))
        if os.path.exists(self.location['original configuration']):
            command = [ 'pheniqs', 'mux', '--compile' ]
            command.append('--config')
            command.append(self.location['original configuration'])

            if self.instruction['sense_input']:
                command.append('--sense-input')

            if 'base_input' in self.instruction:
                command.append('--base-input')
                command.append(self.instruction['base_input'])

            if 'base_output' in self.instruction:
                command.append('--base-output')
                command.append(self.instruction['base_output'])

            if 'input' in self.instruction:
                for value in self.instruction['input']:
                    command.append('--input')
                    command.append(value)

            process = Popen(
                args=command,
                cwd=self.current_working_directoy,
                stdout=PIPE,
                stderr=PIPE
            )
            output, error = process.communicate()
            return_code = process.returncode

            if return_code == 0:
                compiled = json.loads(output.decode('utf8'))
                self.ontology['compiled configurtion'] = compiled
            else:
                raise BadConfigurationError('pheniqs returned {} when creating a static configuration'.format(return_code))
                for line in error.decode('utf8').splitlines():
                    print(error.decode('utf8'))
        else:
            raise BadConfigurationError('could not find configuration file {}'.format(self.instruction['configurtion']))

    def load_report(self):
        # if a report was provided load it
        if 'report' in self.instruction:
            self.location['report'] = os.path.realpath(os.path.abspath(os.path.expanduser(os.path.expandvars(self.instruction['report']))))
            if os.path.exists(self.location['report']):
                with io.open(self.location['report'], 'rb') as file:
                    self.ontology['estimation report'] = json.loads(file.read().decode('utf8'))
            else:
                raise BadConfigurationError('could not find report file {}'.format(self.instruction['report']))
        else:
            self.make_estimation_configuration()
            self.location['estimation configurtion'] = os.path.join(self.current_working_directoy, '{}_estimation_configurtion.json'.format(self.instruction['prefix']))
            with io.open(self.location['estimation configurtion'], 'w') as file:
                file.write(to_json(self.ontology['estimation configurtion']))

            command = [ 'pheniqs', 'mux' ]
            command.append('--config')
            command.append(self.location['estimation configurtion'])

            if self.instruction['sense_input']:
                command.append('--sense-input')

            process = Popen(
                args=command,
                cwd=self.current_working_directoy,
                stdout=PIPE,
                stderr=PIPE
            )
            output, error = process.communicate()
            return_code = process.returncode

            if return_code == 0:
                self.ontology['estimation report'] = json.loads(output.decode('utf8'))
            else:
                print(error)
                raise BadConfigurationError('pheniqs returned {} when estimating prior'.format(return_code))
                for line in error.decode('utf8').splitlines():
                    print(error.decode('utf8'))

    def make_estimation_configuration(self):
        self.ontology['estimation configurtion'] = deepcopy(self.ontology['static configurtion'])
        configuration = self.ontology['estimation configurtion']
        queue = [ configuration ]
        if 'multiplex' in configuration:
            queue += [ configuration['multiplex'] ]
            if 'undetermined' in configuration['multiplex']:
                queue += [ configuration['multiplex']['undetermined'] ]
            if 'codec' in configuration['multiplex']:
                for barcode in configuration['multiplex']['codec'].values():
                    queue += [ barcode ]
        for node in queue:
            if 'output' in node:
                del node['output']

        self.ontology['estimation configurtion']['output'] = [ '/dev/null' ]
        self.ontology['estimation configurtion']['report url'] = '/dev/stdout'

        if 'base_input' in self.instruction:
            self.ontology['estimation configurtion']['base input url'] = self.instruction['base_input']

        if 'base_output' in self.instruction:
            self.ontology['estimation configurtion']['base output url'] = self.instruction['base_output']

        if 'input' in self.instruction:
            self.ontology['estimation configurtion']['input'] = []
            for value in self.instruction['input']:
                self.ontology['estimation configurtion']['input'].append(value)

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    job = None

    try:
        command = CommandLineParser('estimate prior')
        if command.help_triggered:
            command.help()
            sys.exit(0)
        else:
            if 'verbosity' in command.instruction and command.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

            job = EstimatePrior(command.configuration)
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
