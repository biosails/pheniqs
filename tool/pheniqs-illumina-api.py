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
from copy import deepcopy

from core.error import *
from core import log_levels
from core import CommandLineParser
from core import Job
from core import to_json

class IlluminaApi(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)

        if 'preset' in self.instruction and self.instruction['preset'] in self.ontology['preset']:
            self.ontology['selected preset'] = deepcopy(self.ontology['preset'][self.instruction['preset']])
        else:
            self.ontology['selected preset'] = deepcopy(self.ontology['preset']['default'])

    def execute(self):
        self.load_illumina()
        queue = None
        if self.action == 'basecall':
            command = self.write_bcl2fastq_command()
            self.log.info(command)

        elif self.action == 'core':
            # print(to_json(self.instruction))
            # exit(0)
            queue = self.write_core_instruction()

        elif self.action == 'interleave':
            queue = self.write_interleave_instruction_per_lane()

        elif self.action == 'multiplex':
            queue = self.write_multiplex_instruction_per_lane()

        elif self.action == 'estimate':
            queue = self.write_prior_estimate_instruction_per_lane()

        # if queue is not None:
        #     print(to_json(queue))

    @property
    def illumina(self):
        return self.instruction['illumina']

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

    def load_illumina(self):
        self.instruction['illumina'] = {}
        self.parse_run_info()
        self.parse_run_parameters()
        self.parse_sample_sheet()
        self.assemble_platform_model()
        self.assemble_sample_sheet_barcode()
        self.assemble_sample_sheet_lane()

    def parse_run_info(self):
        if 'illumina_run_directory' in self.instruction:
            path = os.path.join(self.instruction['illumina_run_directory'], 'RunInfo.xml')
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
                            self.illumina['DT'] = run_date.isoformat()

                    # flowcell
                    self.illumina['flowcell id'] = run.find('Flowcell').text

                    flowcell_layout = run.find('FlowcellLayout')
                    self.illumina['lane cardinality'] = int(flowcell_layout.attrib['LaneCount'])

                    # instrument
                    self.illumina['instrument id'] = run.find('Instrument').text

                    # read layout
                    self.illumina['segment'] = []
                    index_segment_count = 0
                    not_index_segment_count = 0
                    for read in run.find('Reads'):
                        element = { 'is index': False }
                        is_index = read.attrib['IsIndexedRead']
                        if is_index == 'Y': element['is index'] = True

                        element['cycle cardinality'] = int(read.attrib['NumCycles'])
                        element['illumina segment index'] = int(read.attrib['Number'])
                        element['index'] = element['illumina segment index'] - 1

                        if element['is index']:
                            index_segment_count += 1
                            element['illumina segment name'] = 'I{}'.format(index_segment_count)
                        else:
                            not_index_segment_count += 1
                            element['illumina segment name'] = 'R{}'.format(not_index_segment_count)

                        self.illumina['segment'].append(element)

                    self.illumina['template segment'] = []
                    self.illumina['index segment'] = []
                    for segment in self.illumina['segment']:
                        if segment['is index']:
                            self.illumina['index segment'].append(segment)
                        else:
                            self.illumina['template segment'].append(segment)

                except Exception as error:
                    self.log.warning('failed to parse RunInfo.xml ' + str(error))
            else:
                self.log.warning('RunInfo.xml not found')

    def parse_run_parameters(self):
        if 'illumina_run_directory' in self.instruction:
            path = os.path.join(self.instruction['illumina_run_directory'], 'RunParameters.xml')
            if not os.path.exists(path):
                path = os.path.join(self.instruction['illumina_run_directory'], 'runParameters.xml')
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
                                self.illumina['instrument platform'] = application_name.text.split()[0]

                            application_version = setup.find('ApplicationVersion')
                            if application_version is not None:
                                self.illumina['instrument platform version'] = application_version.text
                except Exception as error:
                    self.log.warning('failed to parse RunParameters.xml ' + str(error))
            else:
                self.log.warning('RunParameters.xml not found')

    def parse_sample_sheet(self):
        if 'illumina_run_directory' in self.instruction:
            path = os.path.join(self.instruction['illumina_run_directory'], 'SampleSheet.csv')
            if os.path.exists(path):
                try:
                    content = None
                    with open(path, 'rb') as file:
                        content = file.read().decode('utf8').splitlines()

                    if content != None:
                        section = None
                        self.illumina['sample sheet'] = {}
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
                                    self.illumina['sample sheet']['header'] = []
                                    pass
                                elif section == 'Reads':
                                    pass
                                elif section == 'Settings':
                                    pass
                                elif section == 'Data':
                                    self.illumina['sample sheet']['data'] = { 'row': [], 'head': [] }

                            else:
                                # process a line in a section
                                if section == 'Header':
                                    line = line.strip()
                                    if line:
                                        self.illumina['sample sheet']['header'].append(line)

                                elif section == 'Reads':
                                    pass
                                elif section == 'Settings':
                                    pass
                                elif section == 'Data':
                                    row = [ x.strip() for x in line.split(',') ]
                                    if len(self.illumina['sample sheet']['data']['head']) == 0:
                                        for key in row:
                                            if key in self.namespace['sample sheet']['Data']['column']:
                                                self.illumina['sample sheet']['data']['head'].append(key)
                                            else:
                                                self.illumina['sample sheet']['data']['head'].append(None)
                                    else:
                                        head = self.illumina['sample sheet']['data']['head']
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
                                            self.illumina['sample sheet']['data']['row'].append(record)
                except Exception as error:
                    self.log.warning('failed to parse SampleSheet.csv ' + str(error))
            else:
                self.log.warning('SampleSheet.csv not found')

    def assemble_platform_model(self):
        PM = None
        if 'instrument platform' in self.illumina:
            PM = self.illumina['instrument platform']

        if 'instrument id' in self.illumina:
            if PM is None:
                PM = self.illumina['instrument id']
            else:
                PM += ' '
                PM += self.illumina['instrument id']

        if PM is not None:
            if 'instrument platform version' in self.illumina:
                PM += ' '
                PM += self.illumina['instrument platform version']

        if PM is not None:
            self.illumina['PM'] = PM

    def assemble_sample_sheet_barcode(self):
        if 'data' in self.illumina['sample sheet'] and self.illumina['sample sheet']['data']['row']:
            for record in self.illumina['sample sheet']['data']['row']:
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

    def assemble_sample_sheet_lane(self):
        if 'data' in self.illumina['sample sheet'] and self.illumina['sample sheet']['data']['row']:
            if all([('lane number' in c) for c in self.illumina['sample sheet']['data']['row']]):
                lane_by_index = {}
                for record in self.illumina['sample sheet']['data']['row']:
                    if record['lane number'] not in lane_by_index:
                        lane_by_index[record['lane number']] = { 'lane number': record['lane number'], 'row': [] }
                    lane_by_index[record['lane number']]['row'].append(record)

                self.illumina['lane'] = []
                for key in sorted(lane_by_index.keys()):
                    self.illumina['lane'].append(lane_by_index[key])

            elif not any([('lane number' in c) for c in self.illumina['sample sheet']['data']['row']]):
                if 'lane cardinality' in self.illumina:
                    lane = { 'lane number': 0, 'row': [] }
                    for record in self.illumina['sample sheet']['data']['row']:
                        lane['row'].append(record)
                    self.illumina['lane'] = [ lane ]
            else:
                raise BadConfigurationError('Incoherent sample sheet, some barcode records define a lane and others dont')

            for lane in self.illumina['lane']:
                # assemble multiplex decoder name
                value = ''
                if 'flowcell id' in self.illumina:
                    value = self.illumina['flowcell id']

                if lane['lane number'] > 0:
                    if len(value) > 0: value += '_'
                    value += 'l{:02d}'.format(lane['lane number'])

                value += '_multiplex'
                lane['multiplex decoder name'] = value

                # assemble multiplex transform
                try:
                    if all(lane['row'][0]['barcode length'] == b['barcode length'] for b in lane['row']):
                        lane['barcode length'] = lane['row'][0]['barcode length']
                        if len(lane['barcode length']) <= len(self.illumina['index segment']):

                            # check index segments have enough cycles for barcode
                            for index,length,segment in zip(range(len(lane['barcode length'])), lane['barcode length'], self.illumina['index segment']):
                                if segment['cycle cardinality'] >= length:
                                    lane['multiplex transform'] = { 'token': [] }
                                    for length,segment in zip(lane['barcode length'], self.illumina['index segment']):
                                        lane['multiplex transform']['token'].append('{}::{}'.format(segment['index'], length))
                                else:
                                    raise BadConfigurationError('read segment ' + segment['index'] + 'has insufficient cycles for barcode segment ' + index)
                        else:
                            raise BadConfigurationError('not enough index read segments for barcode')
                    else:
                        raise BadConfigurationError('not all barcodes are the same length in lane ' + lane['lane number'])

                except BadConfigurationError as error:
                    self.log.warning(str(error))

    def assemble_illumina_templte_transform(self, container):
        if 'template segment' in self.illumina:
            container['transform'] = { 'token': [] }
            for segment in self.illumina['template segment']:
                container['transform']['token'].append('{}::'.format(segment['index']))

    def assemble_illumina_multiplex_decoder(self, lane):
        preset = self.selected_preset['sample sheet record']
        decoder = { 'codec': {} }
        if 'multiplex transform' in lane:
            decoder['transform'] = lane['multiplex transform']

        for record in lane['row']:
            # add inherited global instructions to the barcode record
            for k in [
                'flowcell id',
            ]:
                if k in self.illumina:
                    record[k] = self.illumina[k]

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


    def write_bcl2fastq_command(self):
        bcl2fastq = None
        if 'illumina_run_directory' in self.instruction:
            self.location['basecall samplesheet'] = os.path.join(self.current_working_directoy, 'basecall_samplesheet.csv')
            self.write_basecalling_sample_sheet(self.location['basecall samplesheet'])

            command = [ 'bcl2fastq' ]
            command.append('--runfolder-dir')
            command.append(self.instruction['illumina_run_directory'])
            command.append('--sample-sheet'),
            command.append('basecall_samplesheet.csv')
            command.append('--create-fastq-for-index-reads')
            command.append('--adapter-stringency')
            command.append('0')
            command.append('--minimum-trimmed-read-length')
            command.append('0')
            command.append('--mask-short-adapter-reads')
            command.append('0')

            for key in [
                'no_bgzf_compression',
                'ignore_missing_bcls',
                'ignore_missing_filter',
                'ignore_missing_positions'
            ]:
                if key in self.instruction and self.instruction[key]:
                    command.append('--{}'.format(key.replace('_', '-')))

            if 'output_dir' in self.instruction and self.instruction['output_dir']:
                command.append('--output-dir')
                command.append(self.instruction['output_dir'])

            if 'fastq_compression_level' in self.instruction:
                command.append('--fastq-compression-level')
                command.append(str(self.instruction['fastq_compression_level']))

            command.append('\n')

            bcl2fastq = ' '.join(command)
            path = '{}_basecall.sh'.format(self.illumina['flowcell id'])
            self.log.debug('writing ' + path)
            with io.open(path, 'wb') as file:
                file.write(bcl2fastq.encode('utf8'))

        return bcl2fastq

    def write_basecalling_sample_sheet(self, path):
        buffer = []
        if 'illumina_run_directory' in self.instruction:
            if 'header' in self.illumina['sample sheet']:
                buffer.append('[Header]')
                for line in self.illumina['sample sheet']['header']:
                    buffer.append(line)

            buffer.append('[Data]')
            buffer.append('FCID,Lane,Sample_ID,Sample_Name')
            for lane_number in range(0,self.illumina['lane cardinality']):
                buffer.append('{0},{1},{0},'.format(self.illumina['flowcell id'], lane_number))

        buffer.append('')
        with io.open(path, 'wb') as file:
            file.write('\n'.join(buffer).encode('utf8'))

        self.log.info('writing basecalling sample sheet to ' + path)

    def write_core_instruction(self):
        job = {
            'PL': 'ILLUMINA'
        }
        for key in [
            'DT',
            'PM',
            'flowcell id',
        ]:
            if key in self.illumina:
                job[key] = self.illumina[key]

        for key, name in {
            'base_input_url': 'base input url',
            'base_output_url': 'base output url',
            'no_input_npf': 'filter incoming qc fail',
            'no_output_npf': 'filter outgoing qc fail',
        }.items():
            if key in self.instruction and self.instruction[key]:
                job[name] = self.instruction[key]

        if 'sample sheet' in self.illumina:
            self.assemble_illumina_templte_transform(job)

            if 'lane' in self.illumina and self.illumina['lane']:
                job['decoder'] = {}
                for lane in self.illumina['lane']:
                    job['decoder'][lane['multiplex decoder name']] = self.assemble_illumina_multiplex_decoder(lane)

        path = '{}_core.json'.format(self.illumina['flowcell id'])
        self.log.debug('writing ' + path)
        with io.open(path, 'wb') as file:
            file.write(to_json(job).encode('utf8'))
        return [ job ]

    def write_multiplex_instruction_per_lane(self):
        if 'sample sheet' in self.illumina:
            if 'lane' in self.illumina and self.illumina['lane']:
                queue = []
                core_path = '{}_core.json'.format(self.illumina['flowcell id'])
                for lane in self.illumina['lane']:
                    job = {
                        'import': [ core_path ],
                        'input': [],
                        'output': [],
                        'report url': None,
                        'transform': { 'token': [] },
                        'multiplex': {
                            'algorithm': 'pamld',
                            'base': lane['multiplex decoder name'],
                            'confidence threshold': self.instruction['confidence'],
                            'noise': self.instruction['noise']
                        }
                    }
                    segment_index = 0
                    for segment in self.illumina['segment']:
                        fastq_filename = '{}_S1_L00{}_{}_001.fastq.gz'.format(self.illumina['flowcell id'], lane['lane number'], segment['illumina segment name'])
                        job['input'].append(fastq_filename)
                        if not segment['is index']:
                            job['transform']['token'].append('{}::'.format(segment_index))
                        segment_index += 1

                    job['output'].append('{}_l{:02d}.bam'.format(self.illumina['flowcell id'], lane['lane number']))
                    job['report url'] = '{}_l{:02d}_sample_report.json'.format(self.illumina['flowcell id'], lane['lane number'])

                    path = '{}_l{:02d}_sample.json'.format(self.illumina['flowcell id'], lane['lane number'])
                    self.log.debug('writing ' + path)
                    with io.open(path, 'wb') as file:
                        file.write(to_json(job).encode('utf8'))

                    queue.append(job)
        return queue

    def write_prior_estimate_instruction_per_lane(self):
        if 'sample sheet' in self.illumina:
            if 'lane' in self.illumina and self.illumina['lane']:
                queue = []
                core_path = '{}_core.json'.format(self.illumina['flowcell id'])
                for lane in self.illumina['lane']:
                    job = {
                        'import': [ core_path ],
                        'input': [],
                        'output': [ '/dev/null' ],
                        'report url': None,
                        'transform': { 'token': [] },
                        'multiplex': {
                            'algorithm': 'pamld',
                            'base': lane['multiplex decoder name'],
                            'confidence threshold': self.instruction['confidence'],
                            'noise': self.instruction['noise'],
                            'transform': { 'token': [] }
                        }
                    }
                    segment_index = 0
                    for segment in self.illumina['index segment']:
                        fastq_filename = '{}_S1_L00{}_{}_001.fastq.gz'.format(self.illumina['flowcell id'], lane['lane number'], segment['illumina segment name'])
                        job['input'].append(fastq_filename)
                        token = '{}::'.format(segment_index)
                        job['transform']['token'].append(token)
                        job['multiplex']['transform']['token'].append(token)
                        segment_index += 1

                    job['report url'] = '{}_l{:02d}_estimate_report.json'.format(self.illumina['flowcell id'], lane['lane number'])

                    path = '{}_l{:02d}_estimate.json'.format(self.illumina['flowcell id'], lane['lane number'])
                    self.log.debug('writing ' + path)
                    with io.open(path, 'wb') as file:
                        file.write(to_json(job).encode('utf8'))

                    queue.append(job)
        return queue

    def write_interleave_instruction_per_lane(self):
        if 'sample sheet' in self.illumina:
            if 'lane' in self.illumina and self.illumina['lane']:
                queue = []
                for lane in self.illumina['lane']:
                    job = {
                        'PL': 'ILLUMINA',
                        'input': [],
                        'output': [],
                        'report url': None,
                        'transform': { 'token': [] }
                    }
                    for key in [
                        'DT',
                        'PM',
                        'flowcell id',
                    ]:
                        if key in self.illumina:
                            job[key] = self.illumina[key]

                    segment_index = 0
                    for segment in self.illumina['segment']:
                        fastq_filename = '{}_S1_L00{}_{}_001.fastq.gz'.format(self.illumina['flowcell id'], lane['lane number'], segment['illumina segment name'])
                        job['input'].append(fastq_filename)
                        job['transform']['token'].append('{}::'.format(segment_index))
                        segment_index += 1

                    job['report url'] = '{}_l{:02d}_interleaved_report.json'.format(self.illumina['flowcell id'], lane['lane number'])
                    job['output'].append('{}_l{:02d}_interleaved.bam'.format(self.illumina['flowcell id'], lane['lane number']))

                    path = '{}_l{:02d}_interleaved.json'.format(self.illumina['flowcell id'], lane['lane number'])
                    self.log.debug('writing ' + path)
                    with io.open(path, 'wb') as file:
                        file.write(to_json(job).encode('utf8'))

                    queue.append(job)
        return queue


        job = {
            'PL': 'ILLUMINA'
        }


        if 'segment' in self.illumina:
            for segment in self.illumina['segment']:
                job['transform']['token'].append('{}::'.format(segment['index']))
        self.assemble_interleave_transform(self.interleave)

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    pipeline = None

    try:
        command = CommandLineParser('illumina recipe')
        if command.help_triggered:
            command.help()
            sys.exit(0)
        else:
            if 'verbosity' in command.instruction and command.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

            pipeline = IlluminaApi(command.configuration)
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
