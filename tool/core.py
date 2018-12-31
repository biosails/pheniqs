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
import math
import uuid
import json
import signal
import logging
import hashlib
import platform
from copy import deepcopy
from datetime import datetime, date
from argparse import ArgumentParser
from subprocess import Popen, PIPE

from error import *

log_levels = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

def termination_handler(signal, frame):
    sys.exit(0)

signal.signal(signal.SIGTERM, termination_handler)
signal.signal(signal.SIGINT, termination_handler)

split_class = lambda x: (x[0:x.rfind('.')], x[x.rfind('.') + 1:])
inode_type = lambda path: ( os.path.isfile(path) and 'file' ) or 'directory'
decode_phred = lambda p: ord(p) - 33
encode_phred = lambda q: chr(q + 33)

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

def merge(this, other):
    result = other
    if this is not None:
        result = this
        if other is not None:
            if isinstance(this, dict):
                if isinstance(other, dict):
                    for k,v in other.items():
                        if k in result:
                            result[k] = merge(result[k], v)
                        else:
                            result[k] = v
                else:
                    raise ValueError('incompatible structure')
            else:
                result = other
    return result

def remove_compiled(ontology):
    if ontology is not None:
        if isinstance(ontology, dict):
            for key in list(ontology.keys()):
                if key.endswith('.compiled'):
                    del ontology[key]
                else:
                    next = remove_compiled(ontology[key])
                    if next is None:
                        del ontology[key]
            if len(ontology) == 0:
                ontology = None

        elif isinstance(ontology, list):
            buffer = []
            buffer.extend(ontology)
            ontology.clear()
            for i in buffer:
                next = remove_compiled(i)
                if next is not None:
                    ontology.append(i)

            if len(ontology) == 0:
                ontology = None
    return ontology

def remove_directory(directory, log):
    if os.path.exists(directory):
        log.info('removing {}'.format(directory))
        command = [ 'rm', '-rf' ]
        command.append(directory)
        process = Popen(
            args=command,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate()
        code = process.returncode
        if code is not 0:
            print(output, error, code)
            raise CommandFailedError('failed to remove directory {}'.format(directory))

def prepare_path(path, log, overwrite=True):
    def check_permission(path):
        directory = os.path.dirname(path)
        writable = os.access(directory, os.W_OK)
        present = os.path.exists(directory)

        if writable and present:
            # this hirarchy exists and is writable
            return directory

        elif not (writable or present):
            # try the next one up
            return check_permission(directory)

        elif present and not writable:
            # directory exists but it not writable
            raise PermissionDeniedError(path)

    available = check_permission(path)
    if os.path.exists(path):
        if not overwrite: raise NoOverwriteError(path)
    else:
        directory = os.path.dirname(path)
        if directory != available:
            log.debug('creating directory %s', directory)
            os.makedirs(directory)

def prepare_directory(directory, log):
    def check_permission(directory):
        writable = os.access(directory, os.W_OK)
        present = os.path.exists(directory)

        if writable and present:
            # this hirarchy exists and is writable
            return directory

        elif not (writable or present):
            # try the next one up
            return check_permission(os.path.dirname(directory))

        elif present and not writable:
            # directory exists but it not writable
            raise PermissionDeniedError(directory)

    available = check_permission(directory)
    if available != directory:
       log.debug('creating directory %s', directory)
       os.makedirs(directory)

def f_score(precision, recall):
    if(precision + recall) > 0:
        return 2.0 * (precision * recall) / (precision + recall)
    else:
        return 0

class CommandLineParser(object):
    def __init__(self, name):
        self.ontology = {
            'name': name,
            'configuration path': os.path.join(os.path.dirname(os.path.realpath(__file__)), 'configuration.json'),
            'interface': {},
            'instruction': {}
        }
        self.parser = None
        self.load()

    def load(self):
        with io.open(self.ontology['configuration path'], 'rb') as file:
            content = json.loads(file.read().decode('utf8'))
            self.ontology = merge(self.ontology, content[self.ontology['name']])

        self.parser = ArgumentParser(**self.interface['instruction'])

        # evaluate the type for each prototype
        for prototype in self.interface['prototype'].values():
            if 'type' in prototype['parameter']:
                prototype['parameter']['type'] = eval(prototype['parameter']['type'])

        # add global arguments
        for argument in self.interface['argument']:
            prototype = self.interface['prototype'][argument]
            self.parser.add_argument(*prototype['flag'], **prototype['parameter'])

        if self.sectioned:
            # Add individual command sections
            sub = self.parser.add_subparsers(**self.interface['section']['instruction'])
            for action in self.interface['section']['action']:
                if 'prototype' in action:
                    for prototype in action['prototype'].values():
                        if 'type' in prototype['parameter']:
                            prototype['parameter']['type'] = eval(prototype['parameter']['type'])
                    action['prototype'] = merge(self.interface['prototype'], action['prototype'])
                else:
                    action['prototype'] = deepcopy(self.interface['prototype'])

                key = action['instruction']['name']
                action_parser = sub.add_parser(**action['instruction'])
                if 'argument' in action:
                    for argument in action['argument']:
                        prototype = action['prototype'][argument]
                        action_parser.add_argument(*prototype['flag'], **prototype['parameter'])

                # Add groups of arguments, if any.
                if 'group' in action:
                    for group in action['group']:
                        group_parser = action_parser.add_argument_group(**group['instruction'])
                        if 'argument' in group:
                            for argument in group['argument']:
                                prototype = action['prototype'][argument]
                                group_parser.add_argument(*prototype['flag'], **prototype['parameter'])

        self.ontology['instruction'] = vars(self.parser.parse_args())

    @property
    def help_triggered(self):
        return self.sectioned and self.action is None

    @property
    def interface(self):
        return self.ontology['interface']

    @property
    def sectioned(self):
        return 'section' in self.interface and 'action' in self.interface['section'] and self.interface['section']['action']

    @property
    def instruction(self):
        return self.ontology['instruction']

    @property
    def configuration(self):
        configuration = deepcopy(self.ontology)
        del configuration['interface']
        del configuration['configuration path']
        return configuration

    @property
    def action(self):
        return None if 'action' not in self.instruction else self.instruction['action']

    def help(self):
        self.parser.print_help()

class Job(object):
    def __init__(self, ontology):
        self.log = logging.getLogger('Job')
        self.ontology = {
            'instruction': {
                'home': '~/.pheniqs',
                'platform': platform.system(),
                'current working directoy': os.getcwd(),
            }
        }
        self.ontology = merge(self.ontology, ontology)

    # def __init__(self, name):
    #     self.log = logging.getLogger('Job')
    #     self.ontology = {
    #         'instruction': {
    #             'job name': name,
    #             'platform': platform.system(),
    #             'home': '~/.pheniqs',
    #         }
    #     }
    #
    #     command = CommandLineParser(name)
    #     if command.help_triggered:
    #         command.help()
    #         exit(0)
    #     else:
    #         self.ontology = merge(self.ontology, command.ontology)
    #         del self.ontology['interface']
    #         self.ontology['platform'] = platform.system()
    #
    #         if 'verbosity' in self.instruction and self.instruction['verbosity']:
    #             logging.getLogger().setLevel(log_levels[self.instruction['verbosity']])

    @property
    def instruction(self):
        return self.ontology['instruction']

    @property
    def current_working_directoy(self):
        return self.instruction['current working directoy']

    @property
    def platform(self):
        return self.instruction['platform']

    @property
    def home(self):
        return self.instruction['home']

    @property
    def action(self):
        return self.instruction['action']

    def execute(self):
        pass

    def close(self):
        pass
