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
import sys
import json
import logging
from datetime import datetime

from core.error import *
from core import log_levels
from core import CommandLineParser
from core import Job

class ConfigurationApi(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.banner = [
            'Pheniqs : PHilology ENcoder wIth Quality Statistics',
            'Copyright (C) 2018  Lior Galanti',
            'NYU Center for Genetics and System Biology',
            '',
            'Author: Lior Galanti <lior.galanti@nyu.edu>',
            '',
            'This program is free software: you can redistribute it and/or modify',
            'it under the terms of the GNU Affero General Public License as',
            'published by the Free Software Foundation, either version 3 of the',
            'License, or (at your option) any later version.',
            '',
            'This program is distributed in the hope that it will be useful,',
            'but WITHOUT ANY WARRANTY; without even the implied warranty of',
            'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the',
            'GNU Affero General Public License for more details.',
            '',
            'You should have received a copy of the GNU Affero General Public License',
            'along with this program.  If not, see <http://www.gnu.org/licenses/>.',
        ]
        self.banner.append('')
        self.banner.append('This file is auto generated from {}'.format(self.instruction['path']))
        self.banner.append('')
        self.banner.append('Generated on {}'.format(datetime.now().isoformat()))
        self.configuration = None

    def execute(self):
        if self.action == 'zsh':
            self.execute_zsh_job()

        elif self.action == 'header':
            self.execute_header_job()

    # zsh configuration
    def execute_zsh_job(self):
        if os.path.exists(self.instruction['path']):
            self.log.debug('loading %s', self.instruction['path'])
            with io.open(self.instruction['path'], 'rb') as file:
                self.configuration = json.loads(file.read().decode('utf8'))

            if self.configuration:
                buffer = []
                self.append_zsh_header(buffer)
                self.append_zsh_list_aliases_method(buffer)
                self.append_zsh_commands_method(buffer)

                if 'action' in self.configuration and self.configuration['action']:
                    for action in self.configuration['action']:
                        self.append_zsh_action_method(action, buffer)
                self.append_zsh_main_method(buffer)
                buffer.append('_pheniqs \"$@\"')
                buffer.append('')
                sys.stdout.write('\n'.join(buffer))
        else:
            raise NoConfigurationFileError('could not find configuration {}'.format(self.instruction['path']))

    def append_zsh_header(self, buffer):
        buffer.append('#compdef pheniqs')
        buffer.append('')
        for line in self.banner:
            buffer.append('# {}'.format(line))
        buffer.append('')

    def append_zsh_list_aliases_method(self, buffer):
        buffer.append('_pheniqs_list_aliases() {')
        buffer.append('    local -a aliases')
        buffer.append('    aliases=()')
        buffer.append('    echo "${aliases}"')
        buffer.append('};')
        buffer.append('')

    def append_zsh_commands_method(self, buffer):
        _pheniqs_commands_handler = []
        buffer.append('_pheniqs_commands() {')
        buffer.append('    local -a commands')
        buffer.append('    commands=(')
        if 'action' in self.configuration and self.configuration['action']:
            for action in self.configuration['action']:
                buffer.append('        \'{name}:{description}\''.format(**action))
        buffer.append('    )')
        buffer.append('    _describe -t common-commands \'common commands\' commands')
        buffer.append('};')
        buffer.append('')

    def append_zsh_action_method(self, action, buffer):
        buffer.append('_pheniqs_{name}() {{'.format(**action))
        self.append_zsh_option_arguments(action, buffer)
        buffer.append('};')
        buffer.append('')

    def append_zsh_option_arguments(self, action, buffer):
        buffer.append('    _arguments -C \\')
        if 'option' in action and action['option']:
                for option in action['option']:
                    if self.is_option_optional(option):
                        optspec = [ '    ' ]
                        if self.is_option_plural(option):
                            optspec.append('\*')
                        else:
                            optspec.append('\'({})\''.format(' '.join(option['handle'])))

                        if len(option['handle']) > 1:
                            optspec.append('{{{}}}'.format(','.join(option['handle'])))
                        else:
                            optspec.append('{}'.format(option['handle'][0]))

                        if 'help' in option:
                            optspec.append('\'[{}]\''.format(option['help']))

                        if 'choice' in option:
                            if 'choice description' in option:
                                optspec.append('\':{}:(({}))\''.format(option['name'], ' '.join(['{}\:\"{}\"'.format(c[0], c[1]) for c in zip(option['choice'], option['choice description'])])))
                            else:
                                optspec.append('\':{}:({})\''.format(option['name'], ' '.join(option['choice'])))
                        else:
                            if 'type' in option:
                                if option['type'] == 'boolean':
                                    pass

                                elif option['type'] == 'integer':
                                    optspec.append('\':{name}:\''.format(**option))

                                elif option['type'] == 'decimal':
                                    optspec.append('\':{name}:\''.format(**option))

                                elif option['type'] == 'string':
                                    optspec.append('\':{name}:\''.format(**option))

                                elif option['type'] == 'url':
                                    if 'inode' in option and option['inode'] == 'directory':
                                        optspec.append('\': :_files -/\'')

                                    else:
                                        if 'extension' in option:
                                            if(len(option['extension']) == 1):
                                                optspec.append('\': :_files -g \"*.{}\"\''.format('|'.join(option['extension'])))
                                            else:
                                                optspec.append('\': :_files -g \"*.({})\"\''.format('|'.join(option['extension'])))
                                        else:
                                            optspec.append('\': :\'')

                        optspec.append(' \\')
                        buffer.append(''.join(optspec))

    def is_option_optional(self, option):
        return 'handle' in option and len(option['handle']) > 0

    def is_option_plural(self, option):
        if 'cardinality' in option:
             if isinstance(option['cardinality'], str) and option['cardinality'] == '*':
                 return True

             if isinstance(option['cardinality'], int):
                 if option['cardinality'] > 1:
                     return True
                 else:
                     return False
        else:
            return False

    def append_zsh_main_method(self, buffer):
        buffer.append('_pheniqs(){')
        buffer.append('    local context curcontext="$curcontext" state state_descr line expl')
        buffer.append('    local ret=1')
        self.append_zsh_option_arguments(self.configuration, buffer)
        buffer.append('    \'1:command:->command\' \\')
        buffer.append('    \'*::options:->options\' && return 0')
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

    # configuration header file
    def execute_header_job(self):
        if os.path.exists(self.instruction['path']):
            self.log.debug('loading %s', self.instruction['path'])
            with io.open(self.instruction['path'], 'rb') as file:
                self.configuration = json.loads(file.read().decode('utf8'))

            if self.configuration:
                buffer = []
                self.append_configuration_h_header(buffer)
                buffer.append('')
                buffer.append('#ifndef PHENIQS_CONFIGURATION_H')
                buffer.append('#define PHENIQS_CONFIGURATION_H')
                buffer.append('')

                serial = json.dumps(self.configuration, sort_keys=True, ensure_ascii=True, allow_nan=False,  indent=None)
                binary = [ hex(ord(c)) for c in serial ]
                length = len(binary)
                width = 0xf
                position = 0

                buffer.append('size_t configuration_json_len = {};'.format(length))
                buffer.append('')
                buffer.append('const char configuration_json[] = {')

                while position < length:
                    line = '    {},'.format(', '.join(binary[ position : position + width ]))
                    position += width
                    if position >= length:
                        line = line[:-1]
                    buffer.append(line)

                buffer.append('};')
                buffer.append('')
                buffer.append('#endif /* PHENIQS_CONFIGURATION_H */')
                print('\n'.join(buffer))
        else:
            raise NoConfigurationFileError('could not find configuration {}'.format(self.instruction['path']))

    def append_configuration_h_header(self, buffer):
        buffer.append('/*')
        for line in self.banner:
            buffer.append('    {}'.format(line))
        buffer.append('*/')

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    pipeline = None

    try:
        command = CommandLineParser('configuration api')
        if command.help_triggered:
            command.help()
            sys.exit(0)
        else:
            if 'verbosity' in command.instruction and command.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

            pipeline = ConfigurationApi(command.configuration)
            pipeline.execute()

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
        if pipeline: pipeline.close()

    sys.exit(0)

if __name__ == '__main__':
    main()
