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
import logging
import platform
from copy import deepcopy

from core.error import *
from core import log_levels
from core import merge
from core import remove_compiled

class Job(object):
    def __init__(self, ontology):
        self.log = logging.getLogger('Job')
        default = {
            'instruction': {
                'home': '~/.pheniqs',
                'platform': platform.system(),
                'current working directoy': os.getcwd(),
            },
            'model': {
                'genealogy': {}
            },
            'execution': {}
        }
        self.ontology = merge(default, ontology)
        self.instruction['home'] = os.path.realpath(os.path.abspath(os.path.expanduser(os.path.expandvars(self.instruction['home']))))

        if 'verbosity' in self.instruction and self.instruction['verbosity']:
            self.log.setLevel(log_levels[self.instruction['verbosity']])

    @property
    def home(self):
        return self.instruction['home']

    @property
    def instruction(self):
        return self.ontology['instruction']

    @property
    def model(self):
        return self.ontology['model']

    @property
    def location(self):
        return self.model['location']

    @property
    def genealogy(self):
        return self.model['genealogy']

    @property
    def summary(self):
        if 'summary' not in self.ontology:
            self.ontology['summary'] = deepcopy(self.model)
            remove_compiled(self.ontology['summary'])
        return self.ontology['summary']

    @property
    def execution(self):
        return self.ontology['execution']

    @property
    def current_working_directoy(self):
        return self.instruction['current working directoy']

    @property
    def platform(self):
        return self.instruction['platform']

    @property
    def action(self):
        return self.instruction['action']

    def execute(self):
        pass

    def close(self):
        pass

class Shell(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        default = {
            'execution': {
                'stdout': [],
                'stderr': [],
            }
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
