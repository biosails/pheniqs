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

import urllib.request, urllib.parse, urllib.error
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from http.client import BadStatusLine

from core import *

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
            if 'job implementation' in ontology:
                module, name  = split_class(ontology['job implementation'])
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
                pipeline.log.error('unknown job implementation')

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

class Make(Package):
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

class BZip2(Make):
    def __init__(self, pipeline, node):
        Make.__init__(self, pipeline, node)

    def install_dynamic(self):
        so_basename = 'libbz2.so'
        full_versioned_so_basename = '{}.{}'.format(so_basename, self.version)
        full_versioned_so_package_path = os.path.join(self.package_url, full_versioned_so_basename)
        full_versioned_so_install_path = os.path.join(self.lib_prefix, full_versioned_so_basename)

        self.log.debug('copying %s to %s', full_versioned_so_package_path, full_versioned_so_install_path)
        command = [ 'rsync', full_versioned_so_package_path, full_versioned_so_install_path ]
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

        split_version = self.version.split('.')
        while(len(split_version) > 1):
            split_version = split_version[:-1]
            partial_version = '.'.join(split_version)
            partial_versioned_so_basename = '{}.{}'.format(so_basename, partial_version)
            partial_versioned_so_install_path = os.path.join(self.lib_prefix, partial_versioned_so_basename)

            self.log.info('symlinking %s to %s', full_versioned_so_basename, partial_versioned_so_install_path)
            os.symlink(full_versioned_so_basename, partial_versioned_so_install_path)

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
                        Make.build(self)
                    else:
                        print(code)
                        print(output.decode('utf8'))
                        print(error.decode('utf8'))
                        raise CommandFailedError('make returned {}'.format(code))
        else:
            Make.build(self)

class LibDeflate(Make):
    def __init__(self, pipeline, node):
        Make.__init__(self, pipeline, node)

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

class RapidJSON(Package):
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

class SAMTools(Make):
    def __init__(self, pipeline, node):
        Make.__init__(self, pipeline, node)
        self.node['configure optional'] = [ '--with-htslib={}'.format(self.install_prefix) ]

class PackagePipeline(Pipeline):
    def __init__(self):
        Pipeline.__init__(self, 'package')
        self.package = None
        self.cache = None

    def load_cache(self):
        if 'cache path' in self.instruction:
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
    def package_implementation(self):
        return self.ontology['package implementation']

    @property
    def cache_path(self):
        return self.instruction['cache path']

    @property
    def install_prefix(self):
        return self.instruction['install prefix']

    @property
    def download_prefix(self):
        return self.instruction['download prefix']

    @property
    def package_prefix(self):
        return self.instruction['package prefix']

    @property
    def bin_prefix(self):
        return self.instruction['bin prefix']

    @property
    def include_prefix(self):
        return self.instruction['include prefix']

    @property
    def lib_prefix(self):
        return self.instruction['lib prefix']

    @property
    def filter(self):
        return self.instruction['filter']

    @property
    def force(self):
        return self.instruction['force']

    def execute(self):
        if 'path' in self.instruction:
            self.instruction['path'] = os.path.abspath(os.path.realpath(os.path.expanduser(os.path.expandvars(self.instruction['path']))))
            if os.path.exists(self.instruction['path']):
                self.log.debug('loading %s', self.instruction['path'])
                with io.open(self.instruction['path'], 'rb') as file:
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

                    ontology['document sha1 digest'] = hashlib.sha1(self.instruction['path'].encode('utf8')).hexdigest()
                    self.ontology['instruction'] = merge(self.instruction, ontology)
                    self.load_cache()
                    if self.instruction['document sha1 digest'] not in self.cache['environment']:
                        self.cache['environment'][self.instruction['document sha1 digest']] = { 'package': {} }

                self.persisted_instruction = self.cache['environment'][self.instruction['document sha1 digest']]

                if self.instruction['package']:
                    self.execution['package'] = []
                    prepare_directory(self.home, self.log)
                    prepare_directory(self.install_prefix, self.log)
                    prepare_directory(self.download_prefix, self.log)
                    prepare_directory(self.package_prefix, self.log)
                    self.stdout = io.open(os.path.join(self.home, 'output'), 'a')
                    self.stderr = io.open(os.path.join(self.home, 'error'), 'a')

                    for o in self.instruction['package']:
                        key = o['name']
                        if key in self.package_implementation:
                            o = merge(o, self.package_implementation[key])

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

                                self.save_cache()

    def close(self):
        self.save_cache()
        Pipeline.close(self)

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.ERROR)

    pipeline = None
    try:
        pipeline = PackagePipeline()
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
