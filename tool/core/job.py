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
import json
import logging
import platform
import hashlib
from copy import deepcopy
from datetime import datetime, date

from core.error import *
from core import to_json
from core import log_levels
from core import merge
from core import remove_compiled
from core import prepare_path
from core import prepare_directory

class Job(object):
    def __init__(self, ontology):
        self.log = logging.getLogger('Job')
        default = {
            'instruction': {
                'home': '~/.pheniqs',
                'platform': platform.system(),
                'current working directoy': os.getcwd(),
                'session': None,
            },
            'model': {
                'genealogy': {},
                'location': {},
            },
            'persistence': {
                'model dirty': False,
                'session': None,
                'session sha1': None,
                'voletile session': True,
            },
        }
        self.ontology = merge(default, ontology)
        self.instruction['home'] = os.path.realpath(os.path.abspath(os.path.expanduser(os.path.expandvars(self.instruction['home']))))
        if 'verbosity' in self.instruction and self.instruction['verbosity']:
            self.log.setLevel(log_levels[self.instruction['verbosity']])

    @property
    def instruction(self):
        return self.ontology['instruction']

    @property
    def home(self):
        return self.instruction['home']

    @property
    def action(self):
        return self.instruction['action']

    @property
    def current_working_directoy(self):
        return self.instruction['current working directoy']

    @property
    def platform(self):
        return self.instruction['platform']

    @property
    def model(self):
        return self.ontology['model']

    @property
    def genealogy(self):
        return self.model['genealogy']

    @property
    def location(self):
        return self.model['location']

    @property
    def model_summary(self):
        if 'model summary' not in self.ontology:
            self.ontology['model summary'] = deepcopy(self.model)
            remove_compiled(self.ontology['model summary'])
        return self.ontology['model summary']

    @property
    def session(self):
        if self.ontology['persistence']['session'] is None:
            self.ontology['persistence']['session'] = {
                'created': str(datetime.now())
            }
            self.ontology['persistence']['session sha1']  = hashlib.sha1(to_json(self.ontology['persistence']['session']).encode('utf8')).hexdigest()

            if self.instruction['session'] is not None:
                self.ontology['persistence']['session home'] = os.path.join(self.home, 'session', self.instruction['session'])
                prepare_directory(self.ontology['persistence']['session home'], self.log)

                self.ontology['persistence']['session db path'] = os.path.realpath(os.path.join(self.ontology['persistence']['session home'], 'session.json'))
                prepare_path(self.ontology['persistence']['session db path'], self.log, True)
                self.ontology['persistence']['voletile session'] = False

                if os.path.exists(self.ontology['persistence']['session db path']):
                    self.log.debug('loading session %s', self.instruction['session'])
                    with io.open(path, 'rb') as file:
                        try:
                            content = file.read()
                            self.ontology['persistence']['session'] = json.loads(content.decode('utf8'))
                            self.ontology['persistence']['session sha1'] = hashlib.sha1(content).hexdigest()
                        except json.decoder.JSONDecodeError as e:
                            self.log.warning('ignoring corrupt session %s', self.instruction['session'])
            else:
                self.log.debug('using a voletile session')

        return self.ontology['persistence']['session']

    def save_session(self):
        if not self.ontology['persistence']['voletile session']:
            content = to_json(self.session).encode('utf8')
            checksum  = hashlib.sha1(content).hexdigest()
            if checksum != self.ontology['persistence']['session sha1']:
                # prepare_path(self.ontology['persistence']['session db path'], self.log, True)
                with io.open(self.ontology['persistence']['session db path'], 'wb') as file:
                    self.log.debug('saving session %s', self.instruction['session'])
                    file.write(content)
                self.ontology['persistence']['session sha1'] = checksum
            else:
                self.log.debug('skipping unnecessary session flush')

    @property
    def is_model_dirty(self):
        return self.ontology['persistence']['model dirty']

    @is_model_dirty.setter
    def is_model_dirty(self, value):
        self.ontology['persistence']['model dirty'] = value

    def execute(self):
        pass

    def close(self):
        self.save_session()

class ShellCommand(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.log = logging.getLogger('ShellCommand')
        default = {
            'execution summary': {
                'stdout': [],
                'stderr': [],
            },
            'persistence': {
                'execution summary dirty': False
            },
        }
        self.ontology = merge(default, self.ontology)

        # %E Elapsed real time in [hours:]minutes:seconds
        # %e Elapsed real time in seconds
        # %S Total number of CPU seconds that the process spent in kernel mode
        # %U Total number of CPU seconds that the process spent in user mode
        # %P Percentage of the CPU that this job got computed as (%U + %S) / %E
        # %M Maximum resident set size of the process during its lifetime, in Kbytes
        self.posix_time_command = [ '/usr/local/bin/gtime', '-f', 'POSIX_TIME_START%e,%S,%U,%MPOSIX_TIME_END' ]
        self.posix_time_head_ex = re.compile(r'^POSIX_TIME_START(?P<real>[0-9\.]+),(?P<kernel>[0-9\.]+),(?P<user>[0-9\.]+),(?P<memory>[0-9]+)POSIX_TIME_END$')

    @property
    def execution_summary(self):
        return self.ontology['execution summary']

    @property
    def is_execution_summary_dirty(self):
        return self.ontology['persistence']['execution summary dirty']

    @is_execution_summary_dirty.setter
    def is_execution_summary_dirty(self, value):
        self.ontology['persistence']['execution summary dirty'] = value
