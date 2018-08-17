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

from core import *

class ShellPipeline(Pipeline):
    def __init__(self):
        Pipeline.__init__(self, 'shell')

    def execute(self):
        if self.action == 'zsh':
            self.execute_zsh_job()

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
        if os.path.exists(self.instruction['path']):
            self.log.debug('loading %s', self.instruction['path'])
            buffer = []
            with io.open(self.instruction['path'], 'rb') as file:
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

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)

    pipeline = None
    try:
        pipeline = ShellPipeline()
        pipeline.execute()

    except DownloadError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except ValueError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except CommandFailedError as e:
        logging.getLogger('main').critical(e)
        sys.exit(1)

    except(KeyboardInterrupt, SystemExit) as e:
        if e.code != 0:
            logging.getLogger('main').critical(e)
            sys.exit(1)

    finally:
        if pipeline:
            pipeline.close()

    sys.exit(0)

if __name__ == '__main__':
    main()
