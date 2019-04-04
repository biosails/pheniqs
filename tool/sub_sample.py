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
import numpy
import pickle
import logging
import hashlib
from copy import deepcopy
from datetime import datetime, date
from numpy.random import geometric

from core.error import *
from core import log_levels
from core import CommandLineParser
from core import Job
from core import merge
from core import to_json
from core import prepare_path
from simulation.transcode import Transcode

class SubSampleSam(Transcode):
    def __init__(self, ontology):
        Transcode.__init__(self, ontology)
        self.decay = 0

    def manipulate(self):
        while(len(self.input_buffer) > 0):
            read = self.input_buffer.pop(0)
            if self.decay < 2:
                self.output_buffer.append(read)
                self.decay = geometric(self.instruction['fraction'])
            else:
                self.decay -= 1

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    pipeline = None

    try:
        command = CommandLineParser('subsample')
        if command.help_triggered:
            command.help()
            sys.exit(0)
        else:
            if 'verbosity' in command.instruction and command.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

            pipeline = SubSampleSam(command.configuration)
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
