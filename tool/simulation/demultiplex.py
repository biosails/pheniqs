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
import io
import logging
from subprocess import Popen, PIPE

from core.error import *
from core import ShellCommand
from core import merge
from core import to_json

class PamldDemultiplex(ShellCommand):
    def __init__(self, ontology):
        ShellCommand.__init__(self, ontology)
        self.log = logging.getLogger('PamldDemultiplex')

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    def execute(self):
        self.instruction['input'] = os.path.join(self.home, self.location['simulated substitution path'])
        self.instruction['output'] = os.path.join(self.home, self.location['pamld demultiplex path'])

        if not os.path.exists(self.instruction['output']):
            if os.path.exists(self.instruction['output']):
                self.log.info('purging existing output file %s', self.location['pamld demultiplex path'])
                os.remove(self.instruction['output'])

            command = [ 'pheniqs', 'mux', '--sense-input' ]
            command.append('--config')
            command.append(os.path.join(self.home, self.location['pamld adjusted configuration path']))
            command.append('--input')
            command.append(os.path.join(self.instruction['input']))
            command.append('--output')
            command.append(os.path.join(self.instruction['output']))
            command.append('--report')
            command.append(os.path.join(self.home, self.location['pamld demultiplex report path']))

            self.execution_summary['command'] = ' '.join([str(i) for i in command])
            self.log.debug('executing %s', self.execution_summary['command'])

            process = Popen(
                args=self.posix_time_command + command,
                cwd=self.current_working_directoy,
                stdout=PIPE,
                stderr=PIPE
            )
            output, error = process.communicate()
            self.execution_summary['return code'] = process.returncode

            for line in output.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    self.execution_summary['stdout'].append(line)

            for line in error.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    match = self.posix_time_head_ex.search(line)
                    if match:
                        for k,v in match.groupdict().items():
                            self.execution_summary[k] = float(v)
                    else:
                        self.execution_summary['stderr'].append(line)

            if self.execution_summary['return code'] != 0:
                print(to_json(self.execution_summary))
                raise CommandFailedError('pheniqs returned {} when demultiplexing'.format(self.execution_summary['return code']))
        else:
            self.log.info('skipping pamld demultiplexing because %s exists', self.location['pamld demultiplex path'])

class MddDemultiplex(ShellCommand):
    def __init__(self, ontology):
        ShellCommand.__init__(self, ontology)
        self.log = logging.getLogger('MddDemultiplex')

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    def execute(self):
        self.instruction['input'] = os.path.join(self.home, self.location['simulated substitution path'])
        self.instruction['output'] = os.path.join(self.home, self.location['mdd demultiplex path'])

        if not os.path.exists(self.instruction['output']):
            if os.path.exists(self.instruction['output']):
                self.log.info('purging existing output file %s', self.location['mdd demultiplex path'])
                os.remove(self.instruction['output'])

            command = [ 'pheniqs', 'mux', '--sense-input' ]
            command.append('--config')
            command.append(os.path.join(self.home, self.location['mdd configuration path']))
            command.append('--input')
            command.append(self.instruction['input'])
            command.append('--output')
            command.append(self.instruction['output'])
            command.append('--report')
            command.append(os.path.join(self.home, self.location['mdd demultiplex report path']))

            self.execution_summary['command'] = ' '.join([str(i) for i in command])
            self.log.debug('executing %s', self.execution_summary['command'])

            process = Popen(
                args=self.posix_time_command + command,
                cwd=self.current_working_directoy,
                stdout=PIPE,
                stderr=PIPE
            )
            output, error = process.communicate()
            self.execution_summary['return code'] = process.returncode

            for line in output.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    self.execution_summary['stdout'].append(line)

            for line in error.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    match = self.posix_time_head_ex.search(line)
                    if match:
                        for k,v in match.groupdict().items():
                            self.execution_summary[k] = float(v)
                    else:
                        self.execution_summary['stderr'].append(line)

            if self.execution_summary['return code'] != 0:
                print(to_json(self.execution_summary))
                raise CommandFailedError('pheniqs returned {} when demultiplexing'.format(self.execution_summary['return code']))
        else:
            self.log.info('skipping mdd demultiplexing because %s exists', self.location['mdd demultiplex path'])

class PamldAccuratePriorDemultiplex(ShellCommand):
    def __init__(self, ontology):
        ShellCommand.__init__(self, ontology)
        self.log = logging.getLogger('PamldDemultiplex')

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    def execute(self):
        self.instruction['input'] = os.path.join(self.home, self.location['simulated substitution path'])
        self.instruction['output'] = os.path.join(self.home, self.location['pamld accurate prior demultiplex path'])

        if not os.path.exists(self.instruction['output']):
            if os.path.exists(self.instruction['output']):
                self.log.info('purging existing output file %s', self.instruction['output'])
                os.remove(self.instruction['output'])

            command = [ 'pheniqs', 'mux', '--sense-input' ]
            command.append('--config')
            command.append(os.path.join(self.home, self.location['pamld accurate prior configuration path']))
            command.append('--input')
            command.append(os.path.join(self.instruction['input']))
            command.append('--output')
            command.append(os.path.join(self.instruction['output']))
            command.append('--report')
            command.append(os.path.join(self.home, self.location['pamld accurate prior demultiplex report path']))

            self.execution_summary['command'] = ' '.join([str(i) for i in command])
            self.log.debug('executing %s', self.execution_summary['command'])

            process = Popen(
                args=self.posix_time_command + command,
                cwd=self.current_working_directoy,
                stdout=PIPE,
                stderr=PIPE
            )
            output, error = process.communicate()
            self.execution_summary['return code'] = process.returncode

            for line in output.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    self.execution_summary['stdout'].append(line)

            for line in error.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    match = self.posix_time_head_ex.search(line)
                    if match:
                        for k,v in match.groupdict().items():
                            self.execution_summary[k] = float(v)
                    else:
                        self.execution_summary['stderr'].append(line)

            if self.execution_summary['return code'] != 0:
                print(to_json(self.execution_summary))
                raise CommandFailedError('pheniqs returned {} when demultiplexing'.format(self.execution_summary['return code']))
        else:
            self.log.info('skipping accurate prior pamld demultiplexing because %s exists', self.instruction['output'])

class PamldUniformDemultiplex(ShellCommand):
    def __init__(self, ontology):
        ShellCommand.__init__(self, ontology)
        self.log = logging.getLogger('PamldDemultiplex')

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    def execute(self):
        self.instruction['input'] = os.path.join(self.home, self.location['simulated substitution path'])
        self.instruction['output'] = os.path.join(self.home, self.location['pamld uniform demultiplex path'])

        if not os.path.exists(self.instruction['output']):
            if os.path.exists(self.instruction['output']):
                self.log.info('purging existing output file %s', self.instruction['output'])
                os.remove(self.instruction['output'])

            command = [ 'pheniqs', 'mux', '--sense-input' ]
            command.append('--config')
            command.append(os.path.join(self.home, self.location['pamld uniform configuration path']))
            command.append('--input')
            command.append(os.path.join(self.instruction['input']))
            command.append('--output')
            command.append(os.path.join(self.instruction['output']))
            command.append('--report')
            command.append(os.path.join(self.home, self.location['pamld uniform demultiplex report path']))

            self.execution_summary['command'] = ' '.join([str(i) for i in command])
            self.log.debug('executing %s', self.execution_summary['command'])

            process = Popen(
                args=self.posix_time_command + command,
                cwd=self.current_working_directoy,
                stdout=PIPE,
                stderr=PIPE
            )
            output, error = process.communicate()
            self.execution_summary['return code'] = process.returncode

            for line in output.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    self.execution_summary['stdout'].append(line)

            for line in error.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    match = self.posix_time_head_ex.search(line)
                    if match:
                        for k,v in match.groupdict().items():
                            self.execution_summary[k] = float(v)
                    else:
                        self.execution_summary['stderr'].append(line)

            if self.execution_summary['return code'] != 0:
                print(to_json(self.execution_summary))
                raise CommandFailedError('pheniqs returned {} when demultiplexing'.format(self.execution_summary['return code']))
        else:
            self.log.info('skipping uniform pamld demultiplexing because %s exists', self.instruction['output'])

class DemlDemultiplex(ShellCommand):
    def __init__(self, ontology):
        ShellCommand.__init__(self, ontology)
        self.log = logging.getLogger('DemlDemultiplex')

    @property
    def bsid(self):
        return self.genealogy['bsid']

    @property
    def ssid(self):
        return self.genealogy['ssid']

    def execute(self):
        self.instruction['input'] = os.path.join(self.home, self.location['deml simulated substitution path'])
        self.instruction['output'] = os.path.join(self.home, self.location['deml demultiplex path'])

        if not os.path.exists(self.instruction['output']):
            if os.path.exists(self.instruction['output']):
                self.log.info('purging existing output file %s', self.location['deml demultiplex path'])
                os.remove(self.instruction['output'])

            command = [ 'deML' ]
            command.append('--index')
            command.append(os.path.join(self.home, self.location['deml index path']))
            command.append('--outfile')
            command.append(self.instruction['output'])
            command.append('--summary')
            command.append(os.path.join(self.home, self.location['deml summary path']))
            command.append(self.instruction['input'])

            self.execution_summary['command'] = ' '.join([str(i) for i in command])
            self.log.debug('executing %s', self.execution_summary['command'])

            process = Popen(
                args=self.posix_time_command + command,
                cwd=self.current_working_directoy,
                stdout=PIPE,
                stderr=PIPE
            )
            output, error = process.communicate()
            self.execution_summary['return code'] = process.returncode

            for line in output.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    self.execution_summary['stdout'].append(line)

            for line in error.decode('utf8').splitlines():
                line = line.strip()
                if line:
                    match = self.posix_time_head_ex.search(line)
                    if match:
                        for k,v in match.groupdict().items():
                            self.execution_summary[k] = float(v)
                    else:
                        self.execution_summary['stderr'].append(line)

            if self.execution_summary['return code'] != 0:
                print(to_json(self.execution_summary))
                raise CommandFailedError('deML returned {} when demultiplexing'.format(self.execution_summary['return code']))
        else:
            self.log.info('skipping deML demultiplexing because %s exists', self.location['deml demultiplex path'])
