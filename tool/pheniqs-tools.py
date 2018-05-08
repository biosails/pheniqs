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
from datetime import datetime
from subprocess import Popen, PIPE
from argparse import ArgumentParser
import urllib.request, urllib.parse, urllib.error
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from http.client import BadStatusLine

from error import *
from package import prepare_directory, Package

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

class CommandLineParser(object):
    def __init__(self):
        def add_argument(parser, name):
            node = self.node['prototype'][name]
            parser.add_argument(*node['flag'], **node['parameter'])

        self.ontology = None;
        self.node = None
        configuration_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'configuration.json')
        with io.open(configuration_path, 'rb') as file:
            self.ontology = json.loads(file.read().decode('utf8'))
            self.node = self.ontology['interface']
        self.parser = ArgumentParser(**self.node['instruction'])
        self._instruction = None

        # evaluate the type for each prototype
        for argument in self.node['prototype'].values():
            if 'type' in argument['parameter']:
                argument['parameter']['type'] = eval(argument['parameter']['type'])

        # add global arguments
        for argument in self.node['global']['argument']:
            add_argument(self.parser, argument)

        if self.sectioned:
            # Add individual command sections
            sub = self.parser.add_subparsers(**self.node['section']['instruction'])
            for action in self.node['section']['action']:
                action_parser = sub.add_parser(**action['instruction'])
                if 'argument' in action:
                    for argument in action['argument']:
                        add_argument(action_parser, argument)

                # Add groups of arguments, if any.
                if 'group' in action:
                    for group in action['group']:
                        group_parser = action_parser.add_argument_group(**group['instruction'])
                        if 'argument' in group:
                            for argument in group['argument']:
                                add_argument(group_parser, argument)

    @property
    def sectioned(self):
        return 'section' in self.node and 'action' in self.node['section'] and self.node['section']['action']

    @property
    def instruction(self):
        if self._instruction == None:
            self._instruction = vars(self.parser.parse_args())
        return self._instruction

    @property
    def action(self):
        return None if 'action' not in self.instruction else self.instruction['action']

    def help(self):
        self.parser.print_help()

class Pipeline(object):
    def __init__(self, environment):
        self.log = logging.getLogger('Pipeline')
        self.environment = environment
        self.ontology = self.environment.instruction
        self.execution = {}
        self.stdout = None
        self.stderr = None
        self.package = None
        self.cache = None

    def load_cache(self):
        if 'cache path' in self.ontology:
            if os.path.exists(self.cache_path):
                with io.open(self.cache_path, 'rb') as file:
                    self.cache = json.loads(file.read().decode('utf8'))

            if self.cache is None:
                self.cache = {
                    'environment': {},
                    'created': str(datetime.now()),
                }

            self.cache['loaded'] = str(datetime.now())

    def save_cache(self):
        if 'cache path' in self.ontology:
            self.log.debug('persisting cache')
            with io.open(self.cache_path, 'wb') as file:
                self.cache['saved'] = str(datetime.now())
                content = json.dumps(self.cache, sort_keys=True, ensure_ascii=False, indent=4).encode('utf8')
                file.write(content)

    @property
    def platform(self):
        return self.ontology['platform']

    @property
    def home(self):
        return self.ontology['home']

    @property
    def action(self):
        return self.ontology['action']

    @property
    def cache_path(self):
        return self.ontology['cache path']

    @property
    def install_prefix(self):
        return self.ontology['install prefix']

    @property
    def download_prefix(self):
        return self.ontology['download prefix']

    @property
    def package_prefix(self):
        return self.ontology['package prefix']

    @property
    def bin_prefix(self):
        return self.ontology['bin prefix']

    @property
    def include_prefix(self):
        return self.ontology['include prefix']

    @property
    def lib_prefix(self):
        return self.ontology['lib prefix']

    @property
    def filter(self):
        return self.ontology['filter']

    @property
    def force(self):
        return self.ontology['force']

    def execute(self):
        if self.action in [ 'clean', 'build', 'clean.package' ]:
            self.execute_make_job()

        elif self.action == 'zsh':
            self.execute_zsh_job()

    def close(self):
        self.save_cache()
        if self.stdout:
            self.stdout.close();
        if self.stderr:
            self.stderr.close();

    def execute_make_job(self):
        if 'path' in self.ontology:
            self.ontology['path'] = os.path.abspath(os.path.realpath(os.path.expanduser(os.path.expandvars(self.ontology['path']))))
            if os.path.exists(self.ontology['path']):
                self.log.debug('loading %s', self.ontology['path'])
                with io.open(self.ontology['path'], 'rb') as file:
                    ontology = json.loads(file.read().decode('utf8'))
                    for key in [
                        'home',
                        'platform',
                        'package',
                        'cache path',
                        'install prefix',
                        'download prefix',
                        'package prefix',
                        'bin prefix',
                        'include prefix',
                        'lib prefix',
                    ]:
                        if key not in ontology: ontology[key] = None

                    if not ontology['home']:            ontology['home'] =              '~/.pheniqs'
                    if not ontology['platform']:        ontology['platform'] =          platform.system()
                    if not ontology['cache path']:      ontology['cache path'] =        os.path.join(ontology['home'], 'cache.json')
                    if not ontology['install prefix']:  ontology['install prefix'] =    os.path.join(ontology['home'], 'install')
                    if not ontology['download prefix']: ontology['download prefix'] =   os.path.join(ontology['home'], 'download')
                    if not ontology['package prefix']:  ontology['package prefix'] =    os.path.join(ontology['home'], 'package')
                    if not ontology['bin prefix']:      ontology['bin prefix'] =        os.path.join(ontology['install prefix'], 'bin')
                    if not ontology['include prefix']:  ontology['include prefix'] =    os.path.join(ontology['install prefix'], 'include')
                    if not ontology['lib prefix']:      ontology['lib prefix'] =        os.path.join(ontology['install prefix'], 'lib')

                    for path in [
                        'home',
                        'cache path',
                        'install prefix',
                        'download prefix',
                        'package prefix',
                        'bin prefix',
                        'include prefix',
                        'lib prefix',
                    ]:
                        ontology[path] = os.path.abspath(os.path.expanduser(os.path.expandvars(ontology[path])))

                    ontology['document sha1 digest'] = hashlib.sha1(self.ontology['path'].encode('utf8')).hexdigest()
                    self.ontology = merge(self.environment.instruction, ontology)
                    self.load_cache()
                    if self.ontology['document sha1 digest'] not in self.cache['environment']:
                        self.cache['environment'][self.ontology['document sha1 digest']] = { 'package': {} }

                self.persisted_instruction = self.cache['environment'][self.ontology['document sha1 digest']]

                if self.ontology['package']:
                    self.execution['package'] = []
                    prepare_directory(self.home, self.log)
                    prepare_directory(self.install_prefix, self.log)
                    prepare_directory(self.download_prefix, self.log)
                    prepare_directory(self.package_prefix, self.log)
                    self.stdout = io.open(os.path.join(self.home, 'output'), 'a')
                    self.stderr = io.open(os.path.join(self.home, 'error'), 'a')

                    for o in self.ontology['package']:
                        key = o['name']
                        if self.filter is None or key in self.filter:
                            package = Package.create(self, o)
                            if package:
                                self.execution['package'].append(package)

                                if self.action == 'clean':
                                    self.log.info('cleaning %s', package.display_name)
                                    package.clean()

                                elif self.action == 'build':
                                    if not package.installed:
                                        package.install()
                                    else:
                                        self.log.info('%s is already installed', package.display_name)

                                elif self.action == 'clean.package':
                                    self.log.info('clearing %s', package.display_name)
                                    package.clean_package()

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

class Environment(object):
    def __init__(self):
        self.log = logging.getLogger('Environment')
        self.ontology = {
            'interface': None,
            'instruction': None,
        }

        configuration_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'configuration.json')
        with io.open(configuration_path, 'rb') as file:
            content = json.loads(file.read().decode('utf8'))
            self.ontology = merge(self.ontology, content)

        command = CommandLineParser()
        if command.sectioned and command.action is None:
            command.help()
            exit(0)
        else:
            self.ontology['instruction'] = merge(self.ontology['instruction'], command.instruction)

            if 'verbosity' in self.instruction and self.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[self.instruction['verbosity']])

    @property
    def interface(self):
        return self.ontology['interface']

    @property
    def instruction(self):
        return self.ontology['instruction']

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    environment = Environment()
    pipeline = Pipeline(environment)

    try:
        pipeline.execute()
    except ValueError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)
    except CommandFailedError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)
    except(KeyboardInterrupt, SystemExit) as e:
        pipeline.close()
        sys.exit(1)
    pipeline.close()
    sys.exit(0)

if __name__ == '__main__':
    main()
