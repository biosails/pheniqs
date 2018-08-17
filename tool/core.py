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
import json
import logging
import hashlib
import platform
from copy import deepcopy
from datetime import datetime, date
from argparse import ArgumentParser
from subprocess import Popen, PIPE

from error import *

split_class = lambda x: (x[0:x.rfind('.')], x[x.rfind('.') + 1:])

log_levels = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}
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

            elif isinstance(this, list):
                if isinstance(other, list):
                    result.extend(other)
                else:
                    raise ValueError('incompatible structure')
            else:
                result = other
    return result

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

    def add_argument(self, parser, name):
        prototype = self.interface['prototype'][name]
        parser.add_argument(*prototype['flag'], **prototype['parameter'])

    def load(self):
        with io.open(self.ontology['configuration path'], 'rb') as file:
            content = json.loads(file.read().decode('utf8'))
            self.ontology = merge(self.ontology, content[self.ontology['name']])

        self.parser = ArgumentParser(**self.interface['instruction'])

        # evaluate the type for each prototype
        for argument in self.interface['prototype'].values():
            if 'type' in argument['parameter']:
                argument['parameter']['type'] = eval(argument['parameter']['type'])

        # add global arguments
        for argument in self.interface['global']['argument']:
            self.add_argument(self.parser, argument)

        if self.sectioned:
            # Add individual command sections
            sub = self.parser.add_subparsers(**self.interface['section']['instruction'])
            for action in self.interface['section']['action']:
                key = action['instruction']['name']
                action_parser = sub.add_parser(**action['instruction'])
                if 'argument' in action:
                    for argument in action['argument']:
                        self.add_argument(action_parser, argument)

                # Add groups of arguments, if any.
                if 'group' in action:
                    for group in action['group']:
                        group_parser = action_parser.add_argument_group(**group['instruction'])
                        if 'argument' in group:
                            for argument in group['argument']:
                                self.add_argument(group_parser, argument)

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
    def action(self):
        return None if 'action' not in self.instruction else self.instruction['action']

    def help(self):
        self.parser.print_help()

class Pipeline(object):
    def __init__(self, name):
        self.log = logging.getLogger('Pipeline')
        self.ontology = {
            'instruction': {
                'utility name': name,
                'platform': platform.system(),
                'home': '~/.pheniqs',
            }
        }
        self.execution = {}
        self.stdout = None
        self.stderr = None

        command = CommandLineParser(name)
        if command.help_triggered:
            command.help()
            exit(0)
        else:
            self.ontology = merge(self.ontology, command.ontology)
            del self.ontology['interface']

            self.ontology['platform'] = platform.system()

            if 'verbosity' in self.ontology and self.ontology['verbosity']:
                logging.getLogger().setLevel(log_levels[self.ontology['verbosity']])

    @property
    def instruction(self):
        return self.ontology['instruction']

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
        if self.action == 'zsh':
            self.execute_zsh_job()

    def close(self):
        if self.stdout:
            self.stdout.close();
        if self.stderr:
            self.stderr.close();

    def execute_zsh_job(self):
        def parse_zsh_completion_option(node, buffer):
            if 'option' in node and node['option']:
                    for option in node['option']:
                        argument_handler = []
                        if (len(option['handle']) == 1):
                            argument_handler.append('    \'({}){}[{}]'.format(' '.join(option['handle']),','.join(option['handle']), option['help']))
                        else:
                            argument_handler.append('    \'({})\'{{{}}}\'[{}]'.format(' '.join(option['handle']),','.join(option['handle']), option['help']))
                        if 'choice' in option and option['choice']:
                            if 'choice description' in option and option['choice description']:

                                argument_handler.append(':{}:(({}))'.format(option['name'], ' '.join(['{}\:\"{}\"'.format(c[0], c[1]) for c in zip(option['choice'], option['choice description'])])))
                            else:
                                argument_handler.append(':{}:({})'.format(option['name'], ' '.join(option['choice'])))
                        else:
                            if 'type' in option:
                                if option['type'] == 'boolean':
                                    pass
                                elif option['type'] == 'integer':
                                    argument_handler.append(':{name}:'.format(**option))
                                elif option['type'] == 'decimal':
                                    argument_handler.append(':{name}:'.format(**option))
                                elif option['type'] == 'string':
                                    argument_handler.append(':{name}:'.format(**option))
                                elif option['type'] == 'url':
                                    if 'extension' in option:
                                        if(len(option['extension']) == 1):
                                            argument_handler.append(': :_files -g \"*.{}\"'.format('|'.join(option['extension'])))
                                        else:
                                            argument_handler.append(': :_files -g \"*.({})\"'.format('|'.join(option['extension'])))

                                elif option['type'] == 'directory':
                                        argument_handler.append(': :_files -/')

                        argument_handler.append('\' \\')
                        buffer.append(''.join(argument_handler))
        if os.path.exists(self.ontology['path']):
            self.log.debug('loading %s', self.ontology['path'])
            buffer = []
            with io.open(self.ontology['path'], 'rb') as file:
                node = json.loads(file.read().decode('utf8'))
                buffer.append('#compdef pheniqs')
                buffer.append('')
                buffer.append('# Pheniqs : PHilology ENcoder wIth Quality Statistics')
                buffer.append('# Copyright (C) 2017  Lior Galanti')
                buffer.append('# NYU Center for Genetics and System Biology')
                buffer.append('')
                buffer.append('# Author: Lior Galanti <lior.galanti@nyu.edu>')
                buffer.append('')
                buffer.append('# This program is free software: you can redistribute it and/or modify')
                buffer.append('# it under the terms of the GNU Affero General Public License as')
                buffer.append('# published by the Free Software Foundation, either version 3 of the')
                buffer.append('# License, or (at your option) any later version.')
                buffer.append('')
                buffer.append('# This program is distributed in the hope that it will be useful,')
                buffer.append('# but WITHOUT ANY WARRANTY; without even the implied warranty of')
                buffer.append('# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the')
                buffer.append('# GNU Affero General Public License for more details.')
                buffer.append('')
                buffer.append('# You should have received a copy of the GNU Affero General Public License')
                buffer.append('# along with this program.  If not, see <http://www.gnu.org/licenses/>.')
                buffer.append('')
                buffer.append('_pheniqs_list_aliases() {')
                buffer.append('    local -a aliases')
                buffer.append('    aliases=()')
                buffer.append('    echo "${aliases}"')
                buffer.append('};')
                buffer.append('')

                _pheniqs_commands_handler = []
                _pheniqs_commands_handler.append('_pheniqs_commands() {')
                _pheniqs_commands_handler.append('    local -a commands')
                _pheniqs_commands_handler.append('    commands=(')

                _pheniqs_action_handler = []
                if 'action' in node and node['action']:
                    for action in node['action']:
                        _pheniqs_commands_handler.append('        \'{name}:{description}\''.format(**action))
                        _pheniqs_action_handler.append('_pheniqs_{name}() {{'.format(**action))
                        _pheniqs_action_handler.append('    _arguments \\'.format(**action))
                        parse_zsh_completion_option(action, _pheniqs_action_handler)
                        _pheniqs_action_handler.append('};')
                        _pheniqs_action_handler.append('')

                _pheniqs_commands_handler.append('    )')
                _pheniqs_commands_handler.append('    _describe -t common-commands \'common commands\' commands')
                _pheniqs_commands_handler.append('};')

                buffer.extend(_pheniqs_commands_handler)
                buffer.append('')
                buffer.extend(_pheniqs_action_handler)

                buffer.append('_pheniqs(){')
                buffer.append('    local context curcontext="$curcontext" state state_descr line expl')
                buffer.append('    local ret=1')
                buffer.append('    _arguments -C \\')
                parse_zsh_completion_option(node, buffer)

                buffer.append('        \'1:command:->command\' \\')
                buffer.append('        \'*::options:->options\' && return 0')

                buffer.append('    case \"$state\" in')
                buffer.append('        command) _pheniqs_commands && return 0 ;;')
                buffer.append('        options)')
                buffer.append('            local command_or_alias command')
                buffer.append('            local -A aliases')
                buffer.append('            command_or_alias=\"${line[1]}\"')
                buffer.append('            aliases=($(_pheniqs_list_aliases))')
                buffer.append('            command=\"${aliases[$command_or_alias]:-$command_or_alias}\"')
                buffer.append('            local completion_func=\"_pheniqs_${command//-/_}\"')
                buffer.append('            _call_function ret \"${completion_func}\" && return ret')
                buffer.append('            _message \"a completion function is not defined for command or alias: ${command_or_alias}\"')
                buffer.append('            return 1')
                buffer.append('        ;;')
                buffer.append('    esac')
                buffer.append('};')
                buffer.append('')
                buffer.append('_pheniqs \"$@\"')

            print('\n'.join(buffer))
