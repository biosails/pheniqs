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

class PheniqsIoApi(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.log = logging.getLogger('PheniqsIoApi')
        default = {
            'original': None,
            'static': None,
            'estimating': None,
            'adjusted': None,
            'compiled': None,
        }
        self.ontology = merge(default, self.ontology)
        self.location['original'] = os.path.realpath(os.path.abspath(os.path.expanduser(os.path.expandvars(self.instruction['configuration']))))

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
                    self.ontology['static'] = compiled
                else:
                    raise BadConfigurationError('pheniqs returned {} when creating a static configuration'.format(return_code))
                    for line in error.decode('utf8').splitlines():
                        print(error.decode('utf8'))
            else:
                raise BadConfigurationError('could not find configuration file {}'.format(self.instruction['configurtion']))
        return self.ontology['static']

    @property
    def compiled(self):
        if self.ontology['compiled'] is None:
            if os.path.exists(self.location['original']):
                command = [ 'pheniqs', 'mux', '--compile' ]
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
                    raise BadConfigurationError('pheniqs returned {} when creating a static configuration'.format(return_code))
                    for line in error.decode('utf8').splitlines():
                        print(error.decode('utf8'))
            else:
                raise BadConfigurationError('could not find configuration file {}'.format(self.instruction['configurtion']))
        return self.ontology['compiled']

    @property
    def adjusted(self):
        if self.ontology['adjusted'] is None:
            prefix = self.make_format_prefix()
            suffix = self.make_format_suffix()

            overlay = {}
            if 'multiplex' in self.compiled and self.instruction['split_library']:
                overlay['multiplex'] = {}
                segment_cardinality = self.compiled['multiplex']['segment cardinality']
                if 'undetermined' in self.compiled['multiplex']:
                    overlay['multiplex']['undetermined'] = { 'output': [] }
                    if self.instruction['split_segment'] == 'split':
                        for segment_index in range(1,segment_cardinality + 1):
                            overlay['multiplex']['undetermined']['output'].append('{}_undetermined_s{:0>2}.{}'.format(prefix, segment_index, suffix))
                    else:
                        overlay['multiplex']['undetermined']['output'].append('{}_undetermined.{}'.format(prefix, suffix))

                if 'codec' in self.compiled['multiplex']:
                    overlay['multiplex']['codec'] = {}
                    for barcode_key, barcode in self.compiled['multiplex']['codec'].items():
                        name = self.make_library_name(barcode)
                        overlay['multiplex']['codec'][barcode_key] = { 'output': [] }
                        if self.instruction['split_segment']:
                            for segment_index in range(1,segment_cardinality + 1):
                                overlay['multiplex']['codec'][barcode_key]['output'].append('{}_{}_s{:0>2}.{}'.format(prefix, name, segment_index, suffix))
                        else:
                            overlay['multiplex']['codec'][barcode_key]['output'].append('{}_{}.{}'.format(prefix, name, suffix))
            else:
                overlay['output'] = []
                if self.instruction['split_segment']:
                    segment_cardinality = self.compiled['multiplex']['segment cardinality']
                    for segment_index in range(1,segment_cardinality + 1):
                        overlay['output'].append('{}_s{:0>2}.{}'.format(prefix, segment_index, suffix))
                else:
                    overlay['output'].append('{}.{}'.format(prefix, suffix))

            if self.instruction['static']:
                self.ontology['adjusted'] = deepcopy(self.static)
            else:
                self.ontology['adjusted'] = deepcopy(self.original)

            self.strip_output_directive(self.ontology['adjusted'])
            self.ontology['adjusted'] = merge(self.ontology['adjusted'], overlay)

        return self.ontology['adjusted']

    def execute(self):
        print(to_json(self.adjusted))

    def make_format_suffix(self):
        suffix = self.instruction['format']
        if self.instruction['format'] == 'fastq':
            if 'compression' in self.instruction:
                if self.instruction['compression'] == 'gz' or self.instruction['compression'] == 'bzgf':
                    suffix = '{}.gz'.format(suffix)
            else:
                suffix = '{}.gz'.format(suffix)
        return suffix

    def make_format_prefix(self):
        prefix = None
        if 'prefix' in self.instruction:
            prefix = self.instruction['prefix']
        else:
            if 'flowcell id' in self.compiled:
                prefix = self.compiled['flowcell id']
            else:
                raise BadConfigurationError('must provide prefix if flowcell id is not defined')
        return prefix

    def make_library_name(self, barcode):
        value = None
        if 'LB' in barcode:
            value = barcode['LB']
        elif 'PU' in barcode:
            value = barcode['LB']
        elif 'barcode' in barcode:
            value = ''.join(barcode['barcode'])
        return value

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
        command = CommandLineParser('split output')
        if command.help_triggered:
            command.help()
            sys.exit(0)
        else:
            if 'verbosity' in command.instruction and command.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

            job = PheniqsIoApi(command.configuration)
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
