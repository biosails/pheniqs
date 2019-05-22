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
        self.compile_configuration()
        self.estimate_with_pheniqs()
        self.estimate_prior()

    def remove_compiled_output(self, configuration):
        queue = [ configuration ]

        if 'multiplex' in configuration:
            queue += [ configuration['multiplex'] ]

            if 'undetermined' in configuration['multiplex']:
                queue += [ configuration['multiplex']['undetermined'] ]

            if 'codec' in configuration['multiplex']:
                for barcode in configuration['multiplex']['codec'].values():
                    queue += [ barcode ]

        for node in queue:
            if 'feed' in node:
                del node['feed']

            if 'output' in node:
                del node['output']

    def compile_configuration(self):
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
                self.ontology['compiled'] = compiled
            else:
                raise BadConfigurationError('pheniqs returned {} when compiling configuration'.format(return_code))
                for line in error.decode('utf8').splitlines():
                    print(error.decode('utf8'))
        else:
            raise BadConfigurationError('could not find configuration file {}'.format(self.instruction['configurtion']))

    def estimate_with_pheniqs(self):
        self.ontology['estimation configurtion'] = deepcopy(self.ontology['compiled'])
        self.remove_compiled_output(self.ontology['estimation configurtion'])

        self.location['estimation configurtion'] = os.path.join(self.current_working_directoy, '{}_estimation_configurtion.json'.format(self.instruction['prefix']))
        with io.open(self.location['estimation configurtion'], 'w') as file:
            file.write(to_json(self.ontology['estimation configurtion']))

        command = [ 'pheniqs', 'mux' ]
        command.append('--config')
        command.append(self.location['estimation configurtion'])
        command.append('--output')
        command.append('/dev/null')
        command.append('--report')
        command.append('/dev/stdout')

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
            self.ontology['estimation report'] = json.loads(output.decode('utf8'))
        else:
            raise BadConfigurationError('pheniqs returned {} when estimating prior'.format(return_code))
            for line in error.decode('utf8').splitlines():
                print(error.decode('utf8'))

    def estimate_prior(self):
        def estimate_classifier_prior(classifier_model, classifier_report):
            classifier_model['estimate'] = {}

            # collect featrures
            for item in [
                'count',
                'classified count',
                'classified fraction',
                'average classified confidence',
                'average classified distance',
                'pf count',
                'pf classified count',
                'pf fraction',
                'pf classified fraction',
                'classified pf fraction',
                'average pf classified confidence',
                'average pf classified distance',
                'low confidence count',
                'low conditional confidence count',
            ]:
                if item in classifier_report:
                    classifier_model['estimate'][item] = classifier_report[item]
                else:
                    classifier_model['estimate'][item] = 0

            # estimate noise prior
            noise_count = 0
            if 'low conditional confidence count' in classifier_model['estimate']:
                noise_count += classifier_model['estimate']['low conditional confidence count']

            if 'low confidence count' in classifier_model['estimate']:
                noise_count += (classifier_model['estimate']['low confidence count'] * (1.0 - classifier_model['estimate']['average classified confidence']))

            classifier_model['estimate']['estimated noise'] = noise_count / classifier_model['estimate']['count']
            classifier_model['noise'] = classifier_model['estimate']['estimated noise']
            not_noise = 1.0 - classifier_model['estimate']['estimated noise']
            del classifier_model['estimate']

            barcode_report_by_hash = {}
            for barcode in classifier_report['classified']:
                hash = ''.join(barcode['barcode'])
                barcode_report_by_hash[hash] = barcode

            if 'codec' in classifier_model:
                for barcode_model in classifier_model['codec'].values():
                    hash = ''.join(barcode_model['barcode'])
                    if hash in barcode_report_by_hash:
                        barcode_report = barcode_report_by_hash[hash]
                        barcode_model['estimate'] = {}
                        for item in [
                            'count',
                            'average confidence',
                            'average distance',
                            'pooled classified fraction',
                            'pooled fraction',
                            'pf count',
                            'pf fraction',
                            'average pf confidence',
                            'average pf distance',
                            'pf pooled classified fraction',
                            'pf pooled fraction',
                            'low confidence count',
                        ]:
                            if item in barcode_report:
                                barcode_model['estimate'][item] = barcode_report[item]
                            else:
                                barcode_model['estimate'][item] = 0

                        # barcode prior estimation
                        barcode_model['estimate']['estimated concentration'] = not_noise * barcode_model['estimate']['pf pooled classified fraction']
                        barcode_model['concentration'] = barcode_model['estimate']['estimated concentration']
                        del barcode_model['estimate']

            if 'unclassified' in classifier_model:
                barcode_model = classifier_model['unclassified']
                barcode_report = classifier_report['unclassified']
                barcode_model['estimate'] = {}
                for item in [
                    'count',
                    'pf count',
                    'pf fraction',
                    'pf pooled fraction',
                    'pooled classified fraction',
                    'pooled fraction',
                ]:
                    if item in barcode_report:
                        barcode_model['estimate'][item] = barcode_report[item]
                    else:
                        barcode_model['estimate'][item] = 0

        self.ontology['adjusted configurtion'] = deepcopy(self.ontology['compiled'])
        self.remove_compiled_output(self.ontology['adjusted configurtion'])

        if self.instruction['split_fastq']:
            if 'multiplex' in self.ontology['adjusted configurtion']:
                if 'undetermined' in self.ontology['adjusted configurtion']['multiplex']:
                    barcode = self.ontology['adjusted configurtion']['multiplex']['undetermined']
                    barcode['output'] = []
                    for segment_index in range(1,barcode['segment cardinality'] + 1):
                        barcode['output'].append('{}_undetermined_s{:0>2}.fastq.gz'.format(self.instruction['prefix'], segment_index))

                if 'codec' in self.ontology['adjusted configurtion']['multiplex']:
                    for barcode in self.ontology['adjusted configurtion']['multiplex']['codec'].values():
                        barcode['output'] = []
                        hash = ''.join(barcode['barcode'])
                        for segment_index in range(1,barcode['segment cardinality'] + 1):
                            barcode['output'].append('{}_{}_s{:0>2}.fastq.gz'.format(self.instruction['prefix'], hash, segment_index))

        for classifier_type in [ 'multiplex', 'cellular', 'molecular' ]:
            if classifier_type in self.ontology['estimation report'] and classifier_type in self.ontology['adjusted configurtion']:
                classifier_report = self.ontology['estimation report'][classifier_type]
                classifier_model = self.ontology['adjusted configurtion'][classifier_type]

                if isinstance(classifier_model, dict):
                    estimate_classifier_prior(classifier_model, classifier_report)

                elif isinstance(classifier_model, list):
                    model_by_index = {}
                    for item in classifier_model:
                        model_by_index[item['index']]  = item

                    for report_item in classifier_report:
                        model_item = model_by_index[report_item['index']]
                        estimate_classifier_prior(model_item, report_item)

        self.location['adjusted configurtion'] = os.path.join(self.current_working_directoy, '{}_adjusted.json'.format(self.instruction['prefix']))
        with io.open(self.location['adjusted configurtion'], 'w') as file:
            file.write(to_json(self.ontology['adjusted configurtion']))

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