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

class PermissionDeniedError(Exception):
    def __init__(self, path):
        super(Exception, self).__init__('permission denied for {}'.format(path))
        self.path = path

class NoOverwriteError(Exception):
    def __init__(self, path):
        super(Exception, self).__init__('refusing to overwrite {}'.format(path))
        self.path = path

class DownloadError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__('invalid checksum {}'.format(message))

class CommandFailedError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)

class NoConfigurationFileError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)

class BadConfigurationError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)

class UnsupportedError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)
