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

def estimate_decoder_voodoo_prior(decoder, report):
    if 'average classified confidence' in report:

        # estimate noise prior
        noise_count = 0
        if 'low conditional confidence count' in report:
            noise_count += report['low conditional confidence count']

        if 'low confidence count' in report:
            noise_count += (report['low confidence count'] * (1.0 - report['average classified confidence']))

        decoder['noise'] = noise_count / report['count']
        not_noise = 1.0 - decoder['noise']

        if 'codec' in decoder and 'classified' in report:
            # estimate each barcode prior
            barcode_report_by_hash = {}
            for barcode_report in report['classified']:
                hash = ''.join(barcode_report['barcode'])
                barcode_report_by_hash[hash] = barcode_report

            for barcode_model in decoder['codec'].values():
                hash = ''.join(barcode_model['barcode'])
                if hash in barcode_report_by_hash:
                    barcode_report = barcode_report_by_hash[hash]
                    if 'pf pooled classified fraction' in barcode_report and barcode_report['pf pooled classified fraction'] > 0:
                        barcode_model['concentration'] = not_noise * barcode_report['pf pooled classified fraction']
                    else:
                        barcode_model['concentration'] = 0

class PheniqsPriorApi(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.log = logging.getLogger('PheniqsPriorApi')
        default = {
            'original': None,
            'static': None,
            'estimating': None,
            'adjusted': None,
            'report': None,
        }
        self.ontology = merge(default, self.ontology)
        if 'report' in self.instruction:
            self.location['report'] = os.path.realpath(os.path.abspath(os.path.expanduser(os.path.expandvars(self.instruction['report']))))
        if 'configuration' in self.instruction:
            self.location['original'] = os.path.realpath(os.path.abspath(os.path.expanduser(os.path.expandvars(self.instruction['configuration']))))
        if 'prefix' not in self.instruction and 'flowcell id' in self.static:
            self.instruction['prefix'] = self.static['flowcell id']
        self.location['estimating'] = os.path.join(self.current_working_directoy, '{}_estimation_configurtion.json'.format(self.instruction['prefix']))

    @property
    def original(self):
        if self.ontology['original'] is None:
            if os.path.exists(self.location['original']):
                with io.open(self.location['original'], 'rb') as file:
                    self.ontology['original'] = json.loads(file.read().decode('utf8'))
            else:
                raise BadConfigurationError('could not find configuration file {}'.format(self.instruction['configurtion']))
        return self.ontology['original']

    @property
    def static(self):
        if self.ontology['static'] is None:
            if os.path.exists(self.location['original']):
                command = [ 'pheniqs', 'mux', '--static' ]
                command.append('--config')
                command.append(self.location['original'])

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

                process = Popen (
                    args=command,
                    cwd=self.current_working_directoy,
                    stdout=PIPE,
                    stderr=PIPE
                )
                output, error = process.communicate()
                return_code = process.returncode

                if return_code == 0:
                    self.ontology['static'] = json.loads(output.decode('utf8'))
                else:
                    raise BadConfigurationError('pheniqs returned {} when creating a static configuration'.format(return_code))
                    for line in error.decode('utf8').splitlines():
                        print(error.decode('utf8'))
            else:
                raise BadConfigurationError('could not find configuration file {}'.format(self.instruction['configurtion']))

        return self.ontology['static']

    @property
    def estimating(self):
        if self.ontology['estimating'] is None:
            self.ontology['estimating'] = deepcopy(self.static)
            self.strip_output_directive(self.ontology['estimating'])
            self.ontology['estimating']['output'] = [ '/dev/null' ]
            self.ontology['estimating']['report url'] = '/dev/stdout'
            if 'base_input' in self.instruction:
                self.ontology['estimating']['base input url'] = self.instruction['base_input']
            if 'base_output' in self.instruction:
                self.ontology['estimating']['base output url'] = self.instruction['base_output']
            if 'input' in self.instruction:
                self.ontology['estimating']['input'] = []
                for value in self.instruction['input']:
                    self.ontology['estimating']['input'].append(value)
        return self.ontology['estimating']

    @property
    def adjusted(self):
        if self.ontology['adjusted'] is None:
            self.ontology['adjusted'] = deepcopy(self.static)
            for classifier_type in [ 'multiplex', 'cellular', 'molecular' ]:
                if classifier_type in self.report and classifier_type in self.adjusted:
                    m = self.adjusted[classifier_type]
                    r = self.report[classifier_type]

                    if isinstance(m, dict):
                        estimate_decoder_voodoo_prior(m, r)

                    elif isinstance(m, list):
                        index = 0
                        model_by_index = {}
                        for item in m:
                            item['index'] = index
                            model_by_index[index] = item
                            index += 1

                        for report_item in r:
                            model_item = model_by_index[report_item['index']]
                            estimate_decoder_voodoo_prior(model_item, report_item)

        return self.ontology['adjusted']

    @property
    def report(self):
        if self.ontology['report'] is None:
            if 'report' in self.location:
                if os.path.exists(self.location['report']):
                    with io.open(self.location['report'], 'rb') as file:
                        self.ontology['report'] = json.loads(file.read().decode('utf8'))
                else:
                    raise BadConfigurationError('could not find report file {}'.format(self.instruction['report']))
            else:
                if self.estimating:
                    with io.open(self.location['estimating'], 'w') as file:
                        file.write(to_json(self.estimating))

                    command = [ 'pheniqs', 'mux' ]
                    command.append('--config')
                    command.append(self.location['estimating'])

                    process = Popen(
                        args=command,
                        cwd=self.current_working_directoy,
                        stdout=PIPE,
                        stderr=PIPE
                    )
                    output, error = process.communicate()
                    return_code = process.returncode

                    if return_code == 0:
                        self.ontology['report'] = json.loads(output.decode('utf8'))
                    else:
                        print(error)
                        raise BadConfigurationError('pheniqs returned {} when estimating prior'.format(return_code))
                        for line in error.decode('utf8').splitlines():
                            print(error.decode('utf8'))
        return self.ontology['report']

    def execute(self):
        print(to_json(self.adjusted))

    def strip_output_directive(self, instruction):
        if 'output' in instruction:
            del instruction['output']

        if 'multiplex' in instruction:
            if 'undetermined' in instruction['multiplex']:
                if 'output' in instruction['multiplex']['undetermined']:
                    del instruction['multiplex']['undetermined']['output']

            if 'codec' in instruction['multiplex']:
                for barcode in instruction['multiplex']['codec'].values():
                    if 'output' in barcode:
                        del barcode['output']

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    job = None

    try:
        command = CommandLineParser('prior api')
        if command.help_triggered:
            command.help()
            sys.exit(0)
        else:
            if 'verbosity' in command.instruction and command.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

            job = PheniqsPriorApi(command.configuration)
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
