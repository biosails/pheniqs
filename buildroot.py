#!/usr/bin/env python3

# Pheniqs : PHilology ENcoder wIth Quality Statistics
# Copyright (C) 2017  Lior Galanti
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
from copy import deepcopy
from datetime import datetime
from subprocess import Popen, PIPE
from argparse import ArgumentParser
import urllib.request, urllib.parse, urllib.error
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from http.client import BadStatusLine

base_configuration = {
    'interface': {
        'global': {
            'argument': [
                'version', 
                'verbosity'
            ]
        }, 
        'instruction': {
            'description': 'Lior Galanti lior.galanti@nyu.edu NYU Center for Genomics & Systems Biology'
        }, 
        'prototype': {
            'path': {
                'flag': [
                    'configuration path'
                ], 
                'parameter': {
                    'help': 'configuration path', 
                    'metavar': 'PATH',
                    'nargs': '?'
                }
            }, 
            'filter': {
                'flag': [
                    '-f',
                    '--filter'
                ], 
                'parameter': {
                    'help': 'list of packages', 
                    'metavar': 'PACKAGE',
                    'nargs': '*'
                }
            }, 
            'verbosity': {
                'flag': [
                    '-v', 
                    '--verbosity'
                ], 
                'parameter': {
                    'choices': [
                        'debug', 
                        'info', 
                        'warning', 
                        'error', 
                        'critical'
                    ], 
                    'dest': 'verbosity', 
                    'help': 'logging verbosity level', 
                    'metavar': 'LEVEL'
                }
            }, 
            'version': {
                'flag': [
                    '--version'
                ], 
                'parameter': {
                    'action': 'version', 
                    'version': '%[prog]s 1.0'
                }
            }
        }, 
        'section': {
            'action': [
                {
                    'argument': [
                        'filter',
                        'path'
                    ], 
                    'implementation': 'clean', 
                    'instruction': {
                        'help': 'clean build root environment', 
                        'name': 'clean'
                    }
                },
                {
                    'argument': [
                        'filter',
                        'path'
                    ], 
                    'implementation': 'build', 
                    'instruction': {
                        'help': 'build build root environment', 
                        'name': 'build'
                    }
                }
            ], 
            'instruction': {
                'description': '', 
                'dest': 'action', 
                'help': None, 
                'metavar': 'ACTION', 
                'title': 'pipeline operations'
            }
        }
    }
}
log_levels = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

def to_json(node):
    print(json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4))

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

class PermissionDeniedError(Exception):
    def __init__(self, path):
        super(Exception, self).__init__('permission denied for {}'.format(path))
        self.path = path

class NoOverwriteError(Exception):
    def __init__(self, path):
        super(Exception, self).__init__('refusing to overwrite {}'.format(path))
        self.path = path

class InvalidChecksumError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__('invalid checksum {}'.format(message))

class CommandFailedError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)

class NoConfigurationFileError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)

class CommandLineParser(object):
    def __init__(self, node):
        self.node = node
        self.parser = ArgumentParser(**self.node['instruction'])
        self._instruction = None
        def add_argument(parser, name):
            node = self.node['prototype'][name]
            parser.add_argument(*node['flag'], **node['parameter'])

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
    def __init__(self, ontology):
        self.log = logging.getLogger('Package')
        self.ontology = ontology
        self.stdout = None
        self.stderr = None
        self.package = None
        self.cache = None

        self.load_configuration()
        self.load_cache()
        self.load_packages()
        self.load_work_directory()

    def load_configuration(self):
        if 'configuration path' not in self.ontology or self.ontology['configuration path'] is None:
            self.ontology['configuration path'] = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'buildroot.json')

        if os.path.exists(self.configuration_path):
            self.log.debug('loading %s', self.configuration_path)
            with io.open(self.configuration_path, 'rb') as file:
                node = json.loads(file.read().decode('utf8'))
                for key in [
                    'home',
                    'package',
                    'cache path',
                    'package prefix',
                    'install prefix',
                    'download prefix',
                ]:
                    if key not in node: node[key] = None

                if not node['home']: node['home'] =                         '~/.pheniqs'
                if not node['cache path']: node['cache path'] =             os.path.join(node['home'], 'cache.json')
                if not node['download prefix']: node['download prefix'] =   os.path.join(node['home'], 'download')
                if not node['package prefix']: node['package prefix'] =     os.path.join(node['home'], 'package')
                if not node['install prefix']: node['install prefix'] =     os.path.join(node['home'], 'install')

                for path in [
                    'home',
                    'cache path',
                    'download prefix',
                    'package prefix',
                    'install prefix',
                ]:
                    node[path] = os.path.abspath(os.path.expanduser(os.path.expandvars(node[path])))
                node['configuration digest'] = hashlib.sha1(os.path.abspath(os.path.expanduser(os.path.expandvars(self.configuration_path))).encode('utf8')).hexdigest()

                merge(self.ontology, node)
        else:
            raise NoConfigurationFileError('No buildroot.json configuration found')

    def load_cache(self):
        if os.path.exists(self.cache_path):
            with io.open(self.cache_path, 'rb') as file:
                self.cache = json.loads(file.read().decode('utf8'))

        if self.cache is None:
            self.cache = { 
                'environment': {},
                'created': str(datetime.now()),
            }

        if self.configuration_digest not in self.cache['environment']:
            self.cache['environment'][self.configuration_digest] = {
                'package': {},
            }
        self.cache['loaded'] = str(datetime.now())

    def save_cache(self):
        self.log.debug('persisting cache')
        with io.open(self.cache_path, 'wb') as file:
            self.cache['saved'] = str(datetime.now())
            content = json.dumps(self.cache, sort_keys=True, ensure_ascii=False, indent=4).encode('utf8')
            file.write(content)

    def load_packages(self):
        if self.ontology['package'] is not None:
            self.package = []
            for p in self.ontology['package']:
                key = p['name']
                if self.filter is None or key in self.filter:
                    if   key == 'zlib':
                        package = zlibPackage(self, p)
                    elif key == 'xz':
                        package = xzPackage(self, p)
                    elif key == 'bz2':
                        package = bz2Package(self, p)
                    elif key == 'htslib':
                        package = htslibPackage(self, p)
                    elif key == 'rapidjson':
                        package = rapidjsonPackage(self, p)
                    elif key == 'pheniqs':
                        package = pheniqsPackage(self, p)
                    elif key == 'samtools':
                        package = samtoolsPackage(self, p)
                    else:
                        package = Package(self, p)

                    self.package.append(package)

    def load_work_directory(self):
        prepare_directory(self.home, self.log)
        prepare_directory(self.install_prefix, self.log)
        prepare_directory(self.download_prefix, self.log)
        prepare_directory(self.package_prefix, self.log)
        self.stdout = io.open(os.path.join(self.home, 'output'), 'a')
        self.stderr = io.open(os.path.join(self.home, 'error'), 'a')

    @property
    def configuration_path(self):
        return self.ontology['configuration path']

    @property
    def configuration_digest(self):
        return self.ontology['configuration digest']

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
    def filter(self):
        return self.ontology['filter']

    def execute(self):
        if self.action == 'clean':
            self.clean()

        elif self.action == 'build':
            self.install()

    def close(self):
        self.save_cache()
        self.stdout.close();
        self.stderr.close();

    def install(self):
        for package in self.package:
            if not package.installed:
                self.log.info('installing %s', package.display_name)
                package.install()
            else:
                self.log.info('%s is already installed', package.display_name)

    def clean(self):
        for package in self.package:
            self.log.info('cleaning %s', package.display_name)
            package.clean()

# Package

class Package(object):
    def __init__(self, pipeline, node):
        self.log = logging.getLogger('Package')
        self.pipeline = pipeline
        self.load(node)

    def load(self, node):
        for key in [
            'sha1',
            'version',
            'remote url',
            'package url',
            'download url',
            'extension',
            'compression',
            'display name',
            'path in archive',
            'configuration digest',
        ]:
            if key not in node:
                node[key] = None

        node['display name']
        if node['remote url'] is not None:
            remote_dirname, remote_filename = os.path.split(node['remote url'])
            remote_basename, compression = os.path.splitext(remote_filename)
            remote_basename, extension = os.path.splitext(remote_basename)
            compression = compression.strip('.')
            extension = extension.strip('.')

            if node['compression'] is None:
                node['compression'] = compression

            if node['extension'] is None:
                node['extension'] = extension

            if node['download url'] is None:
                node['download url'] = os.path.join(self.download_prefix, remote_filename)

            if node['path in archive'] is None:
                node['path in archive'] = remote_basename

            if node['version'] is None and remote_basename:
                node['version'] = remote_basename.strip(node['name']).strip('-')

        if node['package url'] is None and node['path in archive'] is not None:
            node['package url'] = os.path.join(self.package_prefix, node['path in archive'])

        if node['display name'] is None:
            node['display name'] = node['name']
            if node['version'] is not None:
                node['display name'] = ' '.join([node['display name'], node['version']])

        for url in [
            'package url',
            'download url'
        ]:
            if node[url] is not None:
                node[url] = os.path.abspath(os.path.expanduser(os.path.expandvars(node[url])))

        content = json.dumps(node, sort_keys=True, ensure_ascii=False)
        node['configuration digest'] = hashlib.sha1(content.encode('utf8')).hexdigest()


        if node['configuration digest'] not in self.cache['package']:
            node['unpacked'] = False
            node['configured'] = False
            node['built'] = False
            node['installed'] = False
            self.cache['package'][node['configuration digest']] = node

        self.node = self.cache['package'][node['configuration digest']]

    @property
    def cache(self):
        return self.pipeline.cache['environment'][self.pipeline.configuration_digest]

    @property
    def name(self):
        return self.node['name']

    @property
    def unpacked(self):
        return self.node['unpacked']

    @property
    def configured(self):
        return self.node['configured']

    @property
    def built(self):
        return self.node['built']

    @property
    def installed(self):
        return self.node['installed']

    @property
    def display_name(self):
        return self.node['display name']

    @property
    def stdout(self):
        return self.pipeline.stdout

    @property
    def stderr(self):
        return self.pipeline.stderr

    @property
    def install_prefix(self):
        return self.pipeline.install_prefix

    @property
    def download_prefix(self):
        return self.pipeline.download_prefix

    @property
    def package_prefix(self):
        return self.pipeline.package_prefix

    @property
    def extension(self):
        return self.node['extension']

    @property
    def compression(self):
        return self.node['compression']

    @property
    def path_in_archive(self):
        return self.node['path in archive']

    @property
    def remote_url(self):
        return self.node['remote url']

    @property
    def download_url(self):
        return self.node['download url']

    @property
    def package_url(self):
        return self.node['package url']

    @property
    def version(self):
        return self.node['version']

    @property
    def sha1(self):
        return self.node['sha1']

    def download(self):
        if self.download_url is not None:
            content = None
            if os.path.exists(self.download_url):
                with open(self.download_url, 'rb') as local:
                    content = local.read()
                    checksum = hashlib.sha1(content).hexdigest()
                    if checksum != self.sha1:
                        self.log.warning('removing corrupt archive %s', self.download_url)
                        os.remove(self.download_url)
                        content = None

            if content is None:
                self.log.debug('fetching %s', self.remote_url)
                request = Request(self.remote_url, None)
                try:
                    response = urlopen(request)
                except BadStatusLine as e:
                    self.log.warning('Bad http status error when requesting %s', self.remote_url)
                except HTTPError as e:
                    self.log.warning('Server returned an error when requesting %s: %s', self.remote_url, e.code)
                except URLError as e:
                    self.log.warning('Could not reach server when requesting %s: %s', self.remote_url, e.reason)
                else:
                    content = response.read()
                    checksum = hashlib.sha1(content).hexdigest()
                    if checksum != self.sha1:
                        raise InvalidChecksumError('{} checksum {} differs from {}'.format(self.display_name, checksum, self.sha1))
                    else:
                        prepare_path(self.download_url, self.log)
                        with open(self.download_url, 'wb') as local:
                            local.write(content)
                        self.log.info('downloaded archive saved %s %s', self.display_name, self.sha1)

    def clean_package(self):
        if self.download_url is not None:
            remove_directory(self.package_url, self.log)
            self.node['unpacked'] = False
            self.node['configured'] = False
            self.node['built'] = False
            self.node['installed'] = False
            # self.pipeline.save_cache()

    def clean(self):
        self.node['configured'] = False
        self.node['built'] = False
        self.node['installed'] = False

    def unpack(self):
        if not self.node['unpacked']:
            self.clean_package()
            if self.download_url is not None:
                self.download()

                self.log.info('unpacking %s', self.display_name)
                command = [ 'tar', '-x' ]
                if self.compression == 'gz':
                    command.append('-z')
                elif self.compression == 'bz2':
                    command.append('-j')
                command.append('-f')
                command.append(self.download_url)
                process = Popen(
                    args=command,
                    cwd=self.package_prefix,
                    stdout=self.stdout,
                    stderr=self.stderr
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['unpacked'] = True
                    # self.pipeline.save_cache()
                else:
                    print(output, error, code)
                    raise CommandFailedError('tar returned {}'.format(code))

    def configure(self):
        if not self.node['configured']:
            self.unpack()
            self.node['configured'] = True

    def build(self):
        if not self.node['built']:
            self.configure()
            self.node['built'] = True

    def install(self):
        if not self.node['installed']:
            self.build()
            self.node['installed'] = True

class makePackage(Package):
    def __init__(self, pipeline, node):
        Package.__init__(self, pipeline, node)

    def clean(self):
        if self.package_url is not None:
            if os.path.exists(os.path.join(self.package_url, 'Makefile')):
                self.log.info('cleaning make environment %s', self.display_name)
                command = [ 'make', 'clean']
                process = Popen(
                    args=command,
                    cwd=self.package_url,
                    stdout=self.stdout,
                    stderr=self.stderr
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['configured'] = False
                    self.node['built'] = False
                    self.node['installed'] = False
                    # self.pipeline.save_cache()
                else:
                    raise CommandFailedError('make clean returned {}'.format(code))
            else:
                self.node['configured'] = False
                self.node['built'] = False
                self.node['installed'] = False

    def configure(self, optional=None):
        if not self.node['configured']:
            if self.package_url is not None:
                self.unpack()
                if os.path.exists(os.path.join(self.package_url, 'configure')):
                    self.log.info('configuring make environment %s', self.display_name)
                    command = [ './configure', '--prefix', self.install_prefix ]
                    if optional: command.extend(optional)
                    process = Popen(
                        args=command,
                        cwd=self.package_url,
                        stdout=self.stdout,
                        stderr=self.stderr
                    )
                    output, error = process.communicate()
                    code = process.returncode
                    if code == 0:
                        self.node['configured'] = True
                        # self.pipeline.save_cache()
                    else:
                        raise CommandFailedError('configure returned {}'.format(code))
                else:
                    self.node['configured'] = True

    def build(self):
        if not self.node['built']:
            self.configure()
            if self.package_url is not None:
                self.log.info('building with make %s', self.display_name)
                command = [ 'make' ]
                process = Popen(
                    args=command,
                    cwd=self.package_url,
                    stdout=self.stdout,
                    stderr=self.stderr
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['built'] = True
                    self.pipeline.save_cache()
                else:
                    raise CommandFailedError('make returned {}'.format(code))

    def install(self):
        if not self.node['installed']:
            self.build()
            if self.package_url is not None:
                self.log.info('installing with make %s', self.display_name)
                command = [ 'make' , 'install' ]
                process = Popen(
                    args=command,
                    cwd=self.package_url,
                    stdout=self.stdout,
                    stderr=self.stderr
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['installed'] = True
                    # self.pipeline.save_cache()
                else:
                    raise CommandFailedError('make install returned {}'.format(code))

class zlibPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

    def clean(self):
        if self.package_url is not None:
            if os.path.exists(os.path.join(self.package_url, 'Makefile')):
                self.log.info('cleaning make environment %s', self.display_name)
                command = [ 'make', 'distclean']
                process = Popen(
                    args=command,
                    cwd=self.package_url,
                    stdout=self.stdout,
                    stderr=self.stderr
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['configured'] = False
                    self.node['built'] = False
                    self.node['installed'] = False
                    # self.pipeline.save_cache()
                else:
                    raise CommandFailedError('make clean returned {}'.format(code))
            else:
                self.node['configured'] = False
                self.node['built'] = False
                self.node['installed'] = False

class xzPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

class bz2Package(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

    def install(self):
        if not self.node['installed']:
            self.build()
            if self.package_url is not None:
                self.log.info('installing with make %s', self.display_name)
                command = [ 'make' , 'install' ]
                command.append('PREFIX={}'.format(self.install_prefix))
                process = Popen(
                    args=command,
                    cwd=self.package_url,
                    stdout=self.stdout,
                    stderr=self.stderr
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['installed'] = True
                else:
                    raise CommandFailedError('make install returned {}'.format(code))

class rapidjsonPackage(Package):
    def __init__(self, pipeline, node):
        Package.__init__(self, pipeline, node)

    def install(self):
        if not self.node['installed']:
            self.build()
            if self.package_url is not None:
                self.log.info('copying header files to include folder %s', self.display_name)
                command = [ 'rsync' , '--recursive' ]
                command.append(os.path.join(self.package_url, 'include/'))
                command.append(os.path.join(self.install_prefix, 'include/'))
                process = Popen(
                    args=command,
                    cwd=self.package_url,
                    stdout=self.stdout,
                    stderr=self.stderr
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['installed'] = True
                else:
                    raise CommandFailedError('make install returned {}'.format(code))

class htslibPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

class samtoolsPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)
    def configure(self, optional=None):
        makePackage.configure(self, [ '--with-htslib={}'.format(self.install_prefix), ])

class pheniqsPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

    def build(self):
        if not self.node['built']:
            self.configure()
            if self.package_url is not None:
                self.log.info('building with make %s', self.display_name)
                command = [ 'make' ]
                command.append('PREFIX={}'.format(self.install_prefix))
                process = Popen(
                    args=command,
                    cwd=self.package_url,
                    stdout=self.stdout,
                    stderr=self.stderr
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['built'] = True
                    self.pipeline.save_cache()
                else:
                    raise CommandFailedError('make returned {}'.format(code))

    def install(self):
        if not self.node['installed']:
            self.build()
            if self.package_url is not None:
                self.log.info('installing with make %s', self.display_name)
                command = [ 'make' , 'install' ]
                command.append('PREFIX={}'.format(self.install_prefix))
                process = Popen(
                    args=command,
                    cwd=self.package_url,
                    stdout=self.stdout,
                    stderr=self.stderr
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['installed'] = True
                else:
                    raise CommandFailedError('make install returned {}'.format(code))

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)

    command = CommandLineParser(base_configuration['interface'])
    if command.sectioned and command.action is None:
        command.help()
    else:
        if 'verbosity' in command.instruction and command.instruction['verbosity']:
            logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

        pipeline = Pipeline(command.instruction)
        try:
            pipeline.execute()
        except ValueError as e:
            logging.getLogger('main').critical(e)
            sys.exit(1)
        except(KeyboardInterrupt, SystemExit) as e:
            pipeline.close()
            sys.exit(1)
        pipeline.close()
    sys.exit(0)

if __name__ == '__main__':
    main()
