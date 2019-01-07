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
import json
import logging
from copy import deepcopy
from argparse import ArgumentParser

from core import merge
from core import to_json

from core.error import *

class CommandLineParser(object):
    def __init__(self, name):
        self.ontology = {
            'name': name,
            'configuration path': os.path.realpath(os.path.join(os.path.dirname(__file__), '../configuration/command.json')),
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

        for k,v in vars(self.parser.parse_args()).items():
            if v is not None:
                self.ontology['instruction'][k] = v

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
