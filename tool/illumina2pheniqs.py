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
from copy import deepcopy

from core.error import *
from core import log_levels
from core import CommandLineParser
from core import Job

class Illumina(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)

        if 'preset' in self.instruction and self.instruction['preset'] in self.ontology['preset']:
            self.ontology['selected preset'] = deepcopy(self.ontology['preset'][self.instruction['preset']])
        else:
            self.ontology['selected preset'] = deepcopy(self.ontology['preset']['default'])

    def execute(self):
        if self.action == 'bcl2fastq':
            self.assemble_bcl2fastq_instruction()
            print(self.bcl2fastq)

        else:
            self.parse_run_info()
            self.parse_run_parameters()
            self.parse_sample_sheet()
            self.assemble_barcode()
            self.assemble_lane()
            self.assemble_lane_multiplex_transform()
            self.assemble_platform_model()

            if self.action == 'core':
                self.assemble_core_instruction()
                print(to_json(self.core))

            if self.action == 'interleave':
                self.assemble_interleave_instruction()
                print(to_json(self.interleave))

            if self.action == 'demultiplex':
                self.assemble_demultiplex_instruction()
                print(to_json(self.demultiplex))

    @property
    def namespace(self):
        return self.ontology['namespace']

    @property
    def core(self):
        return self.ontology['core']

    @property
    def interleave(self):
        return self.ontology['interleave']

    @property
    def demultiplex(self):
        return self.ontology['demultiplex']

    @property
    def bcl2fastq(self):
        return self.ontology['bcl2fastq']

    @property
    def selected_preset(self):
        return self.ontology['selected preset']

    def parse_run_info(self):
        if 'illumina run directory' in self.instruction:
            path = os.path.join(self.instruction['illumina run directory'], 'RunInfo.xml')
            if os.path.exists(path):
                try:
                    import xml.etree.ElementTree
                    tree = xml.etree.ElementTree.parse(path)
                    run = tree.getroot().find('Run')

                    # decode the run date
                    illumina_date_string = run.find('Date').text
                    if illumina_date_string:
                        illumina_date_ex = re.compile(r'^(?P<year>[0-9]{2})(?P<month>[0-9]{2})(?P<day>[0-9]{2})$')
                        match = illumina_date_ex.search(illumina_date_string)
                        if match:
                            illumina_date = dict([(k,int(v)) for k,v in match.groupdict().items()])
                            illumina_date['year'] += 2000
                            run_date = date(**illumina_date)
                            self.instruction['DT'] = run_date.isoformat()

                    # flowcell
                    self.instruction['flowcell id'] = run.find('Flowcell').text

                    flowcell_layout = run.find('FlowcellLayout')
                    self.instruction['lane cardinality'] = int(flowcell_layout.attrib['LaneCount'])

                    # instrument
                    self.instruction['instrument id'] = run.find('Instrument').text

                    # read layout
                    self.instruction['segment'] = []
                    for read in run.find('Reads'):
                        element = { 'is index': False }
                        element['cycle cardinality'] = int(read.attrib['NumCycles'])

                        element['illumina segment index'] = int(read.attrib['Number'])
                        element['index'] = element['illumina segment index'] - 1

                        is_index = read.attrib['IsIndexedRead']
                        if is_index == 'Y': element['is index'] = True

                        self.instruction['segment'].append(element)

                    self.instruction['template segment'] = []
                    self.instruction['index segment'] = []
                    for segment in self.instruction['segment']:
                        if segment['is index']:
                            self.instruction['index segment'].append(segment)
                        else:
                            self.instruction['template segment'].append(segment)

                except Exception as error:
                    self.log.warning('failed to parse RunInfo.xml ' + str(error))
            else:
                self.log.warning('RunInfo.xml not found')

    def parse_run_parameters(self):
        if 'illumina run directory' in self.instruction:
            path = os.path.join(self.instruction['illumina run directory'], 'RunParameters.xml')
            if not os.path.exists(path):
                path = os.path.join(self.instruction['illumina run directory'], 'runParameters.xml')
                if not os.path.exists(path):
                    path = None

            if path is not None:
                try:
                    import xml.etree.ElementTree
                    tree = xml.etree.ElementTree.parse(path)
                    run_parameters = tree.getroot()
                    if run_parameters is not None:
                        setup = run_parameters.find('Setup')
                        if setup is not None:
                            application_name = setup.find('ApplicationName')
                            if application_name is not None:
                                self.instruction['instrument platform'] = application_name.text.split()[0]

                            application_version = setup.find('ApplicationVersion')
                            if application_version is not None:
                                self.instruction['instrument platform version'] = application_version.text
                except Exception as error:
                    self.log.warning('failed to parse RunParameters.xml ' + str(error))
            else:
                self.log.warning('RunParameters.xml not found')

    def parse_sample_sheet(self):
        if 'illumina run directory' in self.instruction:
            path = os.path.join(self.instruction['illumina run directory'], 'SampleSheet.csv')
            if os.path.exists(path):
                try:
                    content = None
                    with open(path, 'rb') as file:
                        content = file.read().decode('utf8').splitlines()

                    if content != None:
                        section = None
                        self.instruction['sample sheet'] = {}
                        section_header_re = re.compile('^\s*\[\s*(?P<section>{})\s*\]'.format('|'.join(self.namespace['sample sheet'].keys())))

                        for line in content:
                            header_match = section_header_re.search(line)
                            if header_match:

                                # wrap up current section
                                if section != None:
                                    if section == 'Header':
                                        pass
                                    elif section == 'Reads':
                                        pass
                                    elif section == 'Settings':
                                        pass
                                    elif section == 'Data':
                                        pass

                                # update current section
                                section = header_match.group('section')

                                # initialize next section
                                if section == 'Header':
                                    pass
                                elif section == 'Reads':
                                    pass
                                elif section == 'Settings':
                                    pass
                                elif section == 'Data':
                                    self.instruction['sample sheet']['data'] = []

                            else:
                                # process a line in a section
                                if section == 'Header':
                                    pass
                                elif section == 'Reads':
                                    pass
                                elif section == 'Settings':
                                    pass
                                elif section == 'Data':
                                    row = [ x.strip() for x in line.split(',') ]
                                    if 'head' not in self.instruction['sample sheet']:
                                        head = []
                                        for key in row:
                                            if key in self.namespace['sample sheet']['Data']['column']:
                                                head.append(key)
                                            else:
                                                head.append(None)
                                        if head:
                                            self.instruction['sample sheet']['head'] = head
                                    else:
                                        head = self.instruction['sample sheet']['head']
                                        record = {}
                                        for index,column in enumerate(row):
                                            if index < len(head) and head[index] != None:
                                                if column:
                                                    record[head[index]] = column

                                        if 'Lane' in record:
                                            try:
                                                record['lane number'] = int(record['Lane'])
                                            except ValueError as e:
                                                pass

                                        if record:
                                            self.instruction['sample sheet']['data'].append(record)
                except Exception as error:
                    self.log.warning('failed to parse SampleSheet.csv ' + str(error))
            else:
                self.log.warning('SampleSheet.csv not found')

    def assemble_barcode(self):
        if 'data' in self.instruction['sample sheet'] and self.instruction['sample sheet']['data']:
            for record in self.instruction['sample sheet']['data']:
                barcode_length = []
                barcode = []
                if 'index' in record and record['index']:
                    barcode.append(record['index'])
                    barcode_length.append(len(record['index']))

                if 'index2' in record and record['index2']:
                    barcode.append(record['index2'])
                    barcode_length.append(len(record['index2']))

                if barcode:
                    record['barcode'] = barcode
                    record['barcode length'] = barcode_length
                    record['concatenated barcode'] = ''.join(barcode)

    def assemble_lane(self):
        if 'data' in self.instruction['sample sheet'] and self.instruction['sample sheet']['data']:
            if all([('lane number' in c) for c in self.instruction['sample sheet']['data']]):
                lane_by_index = {}
                for record in self.instruction['sample sheet']['data']:
                    if record['lane number'] not in lane_by_index:
                        lane_by_index[record['lane number']] = { 'lane number': record['lane number'], 'row': [] }
                    lane_by_index[record['lane number']]['row'].append(record)

                self.instruction['sample sheet']['lane'] = []
                for key in sorted(lane_by_index.keys()):
                    self.instruction['sample sheet']['lane'].append(lane_by_index[key])

            elif not any([('lane number' in c) for c in self.instruction['sample sheet']['data']]):
                if 'lane cardinality' in self.instruction:
                        lane = { 'lane number': 0, 'row': [] }
                        for record in self.instruction['sample sheet']['data']:
                            lane['row'].append(record)
                        self.instruction['sample sheet']['lane'] = [ lane ]
            else:
                raise BadConfigurationError('Incoherent sample sheet, some barcode records define a lane and others dont')

    def assemble_lane_multiplex_transform(self):
        if 'sample sheet' in self.instruction and 'lane' in self.instruction['sample sheet']:
            for lane in self.instruction['sample sheet']['lane']:
                try:
                    if all(lane['row'][0]['barcode length'] == b['barcode length'] for b in lane['row']):
                        lane['barcode length'] = lane['row'][0]['barcode length']
                        if len(lane['barcode length']) <= len(self.instruction['index segment']):

                            # check index segments have enough cycles for barcode
                            for index,length,segment in zip(range(len(lane['barcode length'])), lane['barcode length'], self.instruction['index segment']):
                                if segment['cycle cardinality'] >= length:
                                    lane['multiplex transform'] = { 'token': [] }
                                    for length,segment in zip(lane['barcode length'], self.instruction['index segment']):
                                        lane['multiplex transform']['token'].append('{}::{}'.format(segment['index'], length))
                                else:
                                    raise BadConfigurationError('read segment ' + segment['index'] + 'has insufficient cycles for barcode segment ' + index)
                        else:
                            raise BadConfigurationError('not enough index read segments for barcode')
                    else:
                        raise BadConfigurationError('not all barcodes are the same length in lane ' + lane['lane number'])

                except BadConfigurationError as error:
                    self.log.warning(str(error))

    def assemble_interleave_transform(self, container):
        if 'segment' in self.instruction:
            container['transform'] = {
                "token": []
            }
            for segment in self.instruction['segment']:
                container['transform']['token'].append('{}::'.format(segment['index']))

    def assemble_illumina_demultiplex_transform(self, container):
        if 'template segment' in self.instruction:
            container['transform'] = { "token": [] }
            for segment in self.instruction['template segment']:
                container['transform']['token'].append('{}::'.format(segment['index']))

    def decode_value_by_preset(self, record, preset):
        value = None
        # figure out a key for the record
        if isinstance(preset, str):
            try: value = preset.format(**record)
            except KeyError: pass

        elif isinstance(preset, list):
            for pattern in preset:
                try: value = pattern.format(**record)
                except KeyError: pass
                else: break

        return value

    def assemble_illumina_demultiplex_decoder(self, lane):
        preset = self.selected_preset['sample sheet record']
        decoder = { 'codec': {} }
        if 'multiplex transform' in lane:
            decoder['transform'] = lane['multiplex transform']

        for record in lane['row']:
            # add inherited global instructions to the barcode record
            for k in [
                'flowcell id',
            ]:
                if k in self.instruction:
                    record[k] = self.instruction[k]

            # figure out a key for the record
            key = self.decode_value_by_preset(record, preset['key'])
            if key is not None:
                # if a key was found
                if 'barcode' in record and record['barcode']:
                    # and a barcode sequence
                    element = { 'barcode': record['barcode'] }

                    # decode the rest of the attributes
                    for k,p in preset['value'].items():
                        value = self.decode_value_by_preset(record, p)
                        if value is not None:
                            element[k] = value
                    # and add to the codec
                    decoder['codec'][key] = element

        return decoder

    def assemble_platform_model(self):
        PM = None
        if 'instrument platform' in self.instruction:
            PM = self.instruction['instrument platform']

        if 'instrument id' in self.instruction:
            if PM is None:
                PM = self.instruction['instrument id']
            else:
                PM += ' '
                PM += self.instruction['instrument id']

        if PM is not None:
            if 'instrument platform version' in self.instruction:
                PM += ' '
                PM += self.instruction['instrument platform version']

        if PM is not None:
            self.instruction['PM'] = PM

    def assemble_demultiplex_instruction(self):
        self.ontology['demultiplex'] = {
            'PL': 'ILLUMINA'
        }

        for key in [
            'DT',
            'PM',
            'flowcell id',
        ]:
            if key in self.instruction:
                self.demultiplex[key] = self.instruction[key]

        if 'sample sheet' in self.instruction:
            self.assemble_illumina_demultiplex_transform(self.demultiplex)

            if 'lane' in self.instruction['sample sheet'] and self.instruction['sample sheet']['lane']:
                selected = None
                for lane in self.instruction['sample sheet']['lane']:
                    if lane['lane number'] == self.instruction['lane_number']:
                        selected = lane
                if selected is not None:
                    self.demultiplex['multiplex'] = self.assemble_illumina_demultiplex_decoder(lane)
                else:
                    self.log.warning('lane number %d not found', self.instruction['lane_number'])

    def assemble_core_instruction(self):
        self.ontology['core'] = {
            'PL': 'ILLUMINA'
        }

        for key in [
            'DT',
            'PM',
            'flowcell id',
        ]:
            if key in self.instruction:
                self.core[key] = self.instruction[key]

        if 'sample sheet' in self.instruction:
            self.assemble_illumina_demultiplex_transform(self.core)

            if 'lane' in self.instruction['sample sheet'] and self.instruction['sample sheet']['lane']:
                self.core['decoder'] = {}
                for lane in self.instruction['sample sheet']['lane']:
                    decoder_name = self.infer_decoder_name(lane)
                    self.core['decoder'][decoder_name] = self.assemble_illumina_demultiplex_decoder(lane)

    def assemble_interleave_instruction(self):
        self.ontology['interleave'] = {
            'PL': 'ILLUMINA'
        }

        for key in [
            'DT',
            'PM',
            'flowcell id',
        ]:
            if key in self.instruction:
                self.interleave[key] = self.instruction[key]

        self.assemble_interleave_transform(self.interleave)

    def assemble_bcl2fastq_instruction(self):
        if 'illumina run directory' in self.instruction:
            command = [ 'bcl2fastq' ]
            command.append('-R')
            command.append(self.instruction['illumina run directory'])
            command.append('-o')
            command.append(self.instruction['illumina run directory'])
            command.append('--adapter-stringency')
            command.append('0')
            command.append('--minimum-trimmed-read-length')
            command.append('0')
            command.append('--mask-short-adapter-reads')
            command.append('0')
            command.append('--create-fastq-for-index-reads')
            # command.append('--with-failed-reads')
            self.ontology['bcl2fastq'] = ' '.join(command)

    def infer_decoder_name(self, lane):
        value = ''
        if 'flowcell id' in self.instruction:
            value = self.instruction['flowcell id']

        if lane['lane number'] > 0:
            if len(value) > 0: value += '_'
            value += 'lane_{}'.format(lane['lane number'])

        value += '_multiplex'
        return value

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    pipeline = None

    try:
        command = CommandLineParser('illumina')
        if command.help_triggered:
            command.help()
            sys.exit(0)
        else:
            if 'verbosity' in command.instruction and command.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

            pipeline = Illumina(command.configuration)
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
