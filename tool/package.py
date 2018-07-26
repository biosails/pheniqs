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
from datetime import datetime
from subprocess import Popen, PIPE
import urllib.request, urllib.parse, urllib.error
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from http.client import BadStatusLine

from error import *

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

split_class = lambda x: (x[0:x.rfind('.')], x[x.rfind('.') + 1:])

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

class Package(object):
    def __init__(self, pipeline, node):
        self.log = logging.getLogger('Package')
        self.pipeline = pipeline
        for key in [
            'sha1',
            'version',
            'remote url',
            'package url',
            'download url',
            'remote filename',
            'remote basename',
            'extension',
            'compression',
            'display name',
            'path in archive',
            'document sha1 digest',
        ]:
            if key not in node:
                node[key] = None

        if node['remote url'] is not None:

            if isinstance(node['remote url'], list):
                remote_url = node['remote url'][0]
            else:
                remote_url = node['remote url']

            if node['remote filename'] is None:
                remote_dirname, node['remote filename'] = os.path.split(remote_url)

            remote_basename, compression = os.path.splitext(node['remote filename'])
            if node['compression'] is None and compression:
                node['compression'] = compression.strip('.')

            remote_basename, extension = os.path.splitext(remote_basename)
            if node['extension'] is None and extension:
                node['extension'] = extension.strip('.')

            if node['remote basename'] is None and remote_basename:
                 node['remote basename'] = remote_basename

            if node['download url'] is None:
                node['download url'] = os.path.join(self.download_prefix, node['remote filename'])

            if node['path in archive'] is None:
                node['path in archive'] = node['remote basename']

            if node['version'] is None and node['remote basename']:
                node['version'] = node['remote basename'].strip(node['name']).strip('-')

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
        node['document sha1 digest'] = hashlib.sha1(content.encode('utf8')).hexdigest()

        if node['document sha1 digest'] not in self.pipeline.persisted_instruction['package']:
            node['unpacked'] = False
            node['configured'] = False
            node['built'] = False
            node['installed'] = False
            self.pipeline.persisted_instruction['package'][node['document sha1 digest']] = node

        self.node = self.pipeline.persisted_instruction['package'][node['document sha1 digest']]

    @classmethod
    def create(cls, pipeline, ontology):
        instance = None
        if pipeline and ontology:
            if 'name' in ontology:
                if not('implementation' in ontology and ontology['implementation']):
                    ontology['implementation'] = 'package.{}Package'.format(ontology['name'])

                module, name  = split_class(ontology['implementation'])
                try:
                    implementation_module = __import__(module, fromlist=[name])
                    implementation_class = getattr(implementation_module, name)
                    instance = implementation_class(pipeline, ontology)
                except ImportError as e:
                    pipeline.log.error('no module named %s found when attempting to instantiate %s job implementation', module, ontology['action'])
                    pipeline.log.debug(e)
                except AttributeError as e:
                    pipeline.log.error('class %s not defined in module %s when attempting to instantiate job implementation', name, module)
                    pipeline.log.debug(e)
                except Exception as e:
                    pipeline.log.error('%s %s', type(e), e)
            else:
                pipeline.log.error('make job missing a name')
        return instance

    @property
    def env(self):
        if 'env' not in self.node or self.node['env'] is None:
            self.node['env'] = os.environ.copy()
        return self.node['env']

    @property
    def platform(self):
        return self.pipeline.platform

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
    def bin_prefix(self):
        return self.pipeline.bin_prefix

    @property
    def include_prefix(self):
        return self.pipeline.include_prefix

    @property
    def lib_prefix(self):
        return self.pipeline.lib_prefix

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
            if self.sha1 is None:
                if os.path.exists(self.download_url):
                    self.log.debug('removing old real time archive %s', self.download_url)
                    os.remove(self.download_url)

            if os.path.exists(self.download_url):
                with open(self.download_url, 'rb') as local:
                    content = local.read()
                    checksum = hashlib.sha1(content).hexdigest()
                    if self.sha1 is not None and checksum != self.sha1:
                        self.log.warning('removing corrupt archive %s', self.download_url)
                        os.remove(self.download_url)
                        content = None

            if content is None:
                remote_url = self.remote_url
                if not isinstance(self.remote_url, list):
                    remote_url = [ self.remote_url ]


                error = None
                success = False
                for url in remote_url:
                    self.log.debug('fetching %s', url)
                    request = Request(url, None)
                    try:
                        response = urlopen(request)
                    except BadStatusLine as e:
                        error = 'Bad http status error when requesting {}'.format(url)
                        self.log.warning(error)
                    except HTTPError as e:
                        error = 'Server returned an error when requesting {}: {}'.format(url, e.code)
                        self.log.warning(error)
                    except URLError as e:
                        error = 'Could not reach server when requesting {}: {}'.format(url, e.reason)
                        self.log.warning(error)
                    else:
                        content = response.read()
                        checksum = hashlib.sha1(content).hexdigest()
                        if self.sha1 is not None and checksum != self.sha1:
                            error = '{} checksum {} differs from {}'.format(self.display_name, checksum, self.sha1)
                            self.log.warning(error)
                        else:
                            prepare_path(self.download_url, self.log)
                            with open(self.download_url, 'wb') as local:
                                local.write(content)
                            self.log.info('downloaded archive saved %s %s', self.display_name, self.sha1)
                            success = True
                            break

                if not success:
                    raise DownloadError(error)

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
                command = []
                if self.compression in [ 'gz', 'bz2']:
                    command.append('tar')
                    command.append('-x')
                    if self.compression == 'gz':
                        command.append('-z')
                    elif self.compression == 'bz2':
                        command.append('-j')
                    command.append('-f')
                    command.append(self.download_url)

                elif self.compression in [ 'zip' ]:
                    command.append('unzip')
                    command.append('-x')
                    command.append(self.download_url)

                process = Popen(
                    args=command,
                    env=self.env,
                    cwd=self.package_prefix,
                    stdout=PIPE,
                    stderr=PIPE
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['unpacked'] = True
                    self.stdout.write(output.decode('utf8'))
                    self.stderr.write(error.decode('utf8'))
                    # self.pipeline.save_cache()
                else:
                    print(code)
                    print(output.decode('utf8'))
                    print(error.decode('utf8'))
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
        for key in [
            'configure optional',
            'make build optional',
            'make build target',
        ]:
            if key not in node:
                self.node[key] = None

        if 'make install target' not in self.node:
            self.node['make install target'] = 'install'

        if 'make clean target' not in self.node:
            self.node['make clean target'] = 'clean'

        if 'include prefix in make' not in self.node:
            self.node['include prefix in make'] = False

        self.env['CFLAGS'] = '-I{}'.format(self.include_prefix)
        self.env['LDFLAGS'] = '-L{}'.format(self.lib_prefix)

    @property
    def configure_optional(self):
        return self.node['configure optional']

    @property
    def make_build_optional(self):
        return self.node['make build optional']

    @property
    def make_build_target(self):
        return self.node['make build target']

    @property
    def make_install_target(self):
        return self.node['make install target']

    @property
    def make_clean_target(self):
        return self.node['make clean target']

    @property
    def include_prefix_in_make(self):
        return self.node['include prefix in make']

    def clean(self):
        if self.package_url is not None:
            if os.path.exists(os.path.join(self.package_url, 'Makefile')) and self.make_clean_target:
                self.log.info('cleaning make environment %s', self.display_name)
                command = [ 'make', self.make_clean_target ]

                if self.include_prefix_in_make:
                    command.append('PREFIX={}'.format(self.install_prefix))

                self.log.debug(' '.join([str(i) for i in command]))

                process = Popen(
                    args=command,
                    env=self.env,
                    cwd=self.package_url,
                    stdout=PIPE,
                    stderr=PIPE
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['configured'] = False
                    self.node['built'] = False
                    self.node['installed'] = False
                    self.stdout.write(output.decode('utf8'))
                    self.stderr.write(error.decode('utf8'))
                else:
                    print(code)
                    print(output.decode('utf8'))
                    print(error.decode('utf8'))
                    raise CommandFailedError('make clean returned {}'.format(code))
            else:
                self.node['configured'] = False
                self.node['built'] = False
                self.node['installed'] = False

    def configure(self):
        if not self.node['configured']:
            if self.package_url is not None:
                self.unpack()
                if os.path.exists(os.path.join(self.package_url, 'configure')):
                    self.log.info('configuring make environment %s', self.display_name)
                    command = [ './configure', '--prefix={}'.format(self.install_prefix) ]

                    if self.configure_optional:
                        command.extend(self.configure_optional)

                    self.log.debug(' '.join([str(i) for i in command]))

                    process = Popen(
                        args=command,
                        env=self.env,
                        cwd=self.package_url,
                        stdout=PIPE,
                        stderr=PIPE
                    )
                    output, error = process.communicate()
                    code = process.returncode
                    if code == 0:
                        self.node['configured'] = True
                        self.stdout.write(output.decode('utf8'))
                        self.stderr.write(error.decode('utf8'))
                    else:
                        print(code)
                        print(output.decode('utf8'))
                        print(error.decode('utf8'))
                        raise CommandFailedError('configure returned {}'.format(code))
                else:
                    self.node['configured'] = True

    def build(self):
        if not self.node['built']:
            self.configure()
            if self.package_url is not None:
                self.log.info('building with make %s', self.display_name)
                command = [ 'make' ]

                if self.make_build_target:
                    command.append(self.make_build_target)

                if self.include_prefix_in_make:
                    command.append('PREFIX={}'.format(self.install_prefix))

                if self.make_build_optional:
                    command.extend(self.make_build_optional)

                self.log.debug(' '.join([str(i) for i in command]))

                process = Popen(
                    args=command,
                    env=self.env,
                    cwd=self.package_url,
                    stdout=PIPE,
                    stderr=PIPE
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['built'] = True
                    # self.pipeline.save_cache()
                    self.stdout.write(output.decode('utf8'))
                    self.stderr.write(error.decode('utf8'))
                else:
                    print(code)
                    print(output.decode('utf8'))
                    print(error.decode('utf8'))
                    raise CommandFailedError('make returned {}'.format(code))

    def install(self):
        if not self.node['installed']:
            self.build()
            if self.package_url is not None:
                self.log.info('installing with make %s', self.display_name)
                command = [ 'make' , self.make_install_target ]

                if self.include_prefix_in_make:
                    command.append('PREFIX={}'.format(self.install_prefix))

                self.log.debug(' '.join([str(i) for i in command]))

                process = Popen(
                    args=command,
                    env=self.env,
                    cwd=self.package_url,
                    stdout=PIPE,
                    stderr=PIPE
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['installed'] = True
                    self.stdout.write(output.decode('utf8'))
                    self.stderr.write(error.decode('utf8'))
                else:
                    print(code)
                    print(output.decode('utf8'))
                    print(error.decode('utf8'))
                    raise CommandFailedError('make install returned {}'.format(code))

class zlibPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

class xzPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

class bz2Package(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

    def install_dynamic(self):
        dynamic_library_path = os.path.join(self.package_url, 'libbz2.so')
        versioned_dynamic_library_path = '{}.{}'.format(dynamic_library_path, self.version)

        self.log.debug('copying %s to %s', versioned_dynamic_library_path, self.lib_prefix)
        command = [ 'rsync', versioned_dynamic_library_path, self.lib_prefix ]
        process = Popen(
            args=command,
            env=self.env,
            cwd=self.package_url,
            stdout=PIPE,
            stderr=PIPE
        )
        output, error = process.communicate()
        code = process.returncode
        if code == 0:
            self.stdout.write(output.decode('utf8'))
            self.stderr.write(error.decode('utf8'))
        else:
            print(code)
            print(output.decode('utf8'))
            print(error.decode('utf8'))
            raise CommandFailedError('rsync returned {}'.format(code))

        self.log.info('symlinking %s to %s', versioned_dynamic_library_path, dynamic_library_path)
        os.symlink(versioned_dynamic_library_path, dynamic_library_path)

    def build(self):
        if self.platform == 'Linux':
            if not self.node['built']:
                self.configure()
                if self.package_url is not None:
                    self.log.info('building %s dynamic library', self.display_name)
                    command = [ 'make', '--file', 'Makefile-libbz2_so' ]

                    if self.include_prefix_in_make:
                        command.append('PREFIX={}'.format(self.install_prefix))

                    self.log.debug(' '.join([str(i) for i in command]))

                    process = Popen(
                        args=command,
                        env=self.env,
                        cwd=self.package_url,
                        stdout=PIPE,
                        stderr=PIPE
                    )
                    output, error = process.communicate()
                    code = process.returncode
                    if code == 0:
                        self.stdout.write(output.decode('utf8'))
                        self.stderr.write(error.decode('utf8'))
                        self.install_dynamic()
                        makePackage.build(self)
                    else:
                        print(code)
                        print(output.decode('utf8'))
                        print(error.decode('utf8'))
                        raise CommandFailedError('make returned {}'.format(code))
        else:
            makePackage.build(self)

class libdeflatePackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

    def install(self):
        if not self.node['installed']:
            self.build()
            if self.package_url is not None:
                static_library_path = os.path.join(self.package_url, 'libdeflate.a')
                self.log.debug('copying %s to %s', static_library_path, self.lib_prefix)
                command = [ 'rsync', static_library_path, self.lib_prefix ]
                process = Popen(
                    args=command,
                    env=self.env,
                    cwd=self.package_url,
                    stdout=PIPE,
                    stderr=PIPE
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.stdout.write(output.decode('utf8'))
                    self.stderr.write(error.decode('utf8'))
                else:
                    print(code)
                    print(output.decode('utf8'))
                    print(error.decode('utf8'))
                    raise CommandFailedError('rsync returned {}'.format(code))

                library_header_path = os.path.join(self.package_url, 'libdeflate.h')
                self.log.debug('copying %s to %s', library_header_path, self.include_prefix)
                command = [ 'rsync', library_header_path, self.include_prefix ]
                process = Popen(
                    args=command,
                    env=self.env,
                    cwd=self.package_url,
                    stdout=PIPE,
                    stderr=PIPE
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.stdout.write(output.decode('utf8'))
                    self.stderr.write(error.decode('utf8'))
                else:
                    print(code)
                    print(output.decode('utf8'))
                    print(error.decode('utf8'))
                    raise CommandFailedError('rsync returned {}'.format(code))

                if self.platform == 'Linux':
                    dynamic_library_path = os.path.join(self.package_url, 'libdeflate.so')
                    self.log.debug('copying %s to %s', dynamic_library_path, self.lib_prefix)
                    command = [ 'rsync', dynamic_library_path, self.lib_prefix ]
                    process = Popen(
                        args=command,
                        env=self.env,
                        cwd=self.package_url,
                        stdout=PIPE,
                        stderr=PIPE
                    )
                    output, error = process.communicate()
                    code = process.returncode
                    if code == 0:
                        self.stdout.write(output.decode('utf8'))
                        self.stderr.write(error.decode('utf8'))
                    else:
                        print(code)
                        print(output.decode('utf8'))
                        print(error.decode('utf8'))
                        raise CommandFailedError('rsync returned {}'.format(code))

                self.node['installed'] = True

class rapidjsonPackage(Package):
    def __init__(self, pipeline, node):
        Package.__init__(self, pipeline, node)

    def install(self):
        if not self.node['installed']:
            self.build()
            if self.package_url is not None:
                self.log.debug('copying %s header files to %s', self.display_name, self.include_prefix)
                command = [ 'rsync' , '--recursive', os.path.join(self.package_url, 'include/'), self.include_prefix ]
                process = Popen(
                    args=command,
                    env=self.env,
                    cwd=self.package_url,
                    stdout=PIPE,
                    stderr=PIPE
                )
                output, error = process.communicate()
                code = process.returncode
                if code == 0:
                    self.node['installed'] = True
                    self.stdout.write(output.decode('utf8'))
                    self.stderr.write(error.decode('utf8'))
                else:
                    print(code)
                    print(output.decode('utf8'))
                    print(error.decode('utf8'))
                    raise CommandFailedError('rsync returned {}'.format(code))

class htslibPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)

class samtoolsPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)
        self.node['configure optional'] = [ '--with-htslib={}'.format(self.install_prefix) ]

class pheniqsPackage(makePackage):
    def __init__(self, pipeline, node):
        makePackage.__init__(self, pipeline, node)
