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
import sys
import json
import signal
import logging
from copy import deepcopy

from core.error import *
from subprocess import Popen, PIPE

log_levels = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

def to_json(ontology):
    return json.dumps(ontology, sort_keys=True, ensure_ascii=False, indent=4)

def merge(this, other):
    if this is None:
        return deepcopy(other)
    else:
        # this is not None
        if other is None:
            return deepcopy(this)
        else:
            # other is not None
            if isinstance(this, dict):
                if isinstance(other, dict):
                    merged = dict()
                    for k,v in this.items():
                        if k not in other:
                            merged[k] = deepcopy(this[k])
                    for k,v in other.items():
                        if k not in this:
                            merged[k] = deepcopy(v)
                        else:
                            merged[k] = merge(this[k], v)
                    return merged
                else:
                    raise ValueError('incompatible structure')
            else:
                return deepcopy(other)

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
        if code != 0:
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

def split_class(name):
    return (name[0:name.rfind('.')], name[name.rfind('.') + 1:])

def inode_type(path):
    if os.path.isfile(path):
        return 'file'
    else:
        return 'directory'

from core.command_line import CommandLineParser
from core.job import Job

def termination_handler(signal, frame):
    sys.exit(0)

signal.signal(signal.SIGTERM, termination_handler)
signal.signal(signal.SIGINT, termination_handler)
