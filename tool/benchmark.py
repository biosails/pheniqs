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

    # SIMULATION_PRESET="A5KVK"
    #
    # PHENIQS_CODE_HOME="~/code/biosails/pheniqs"
    # DEML_CMD="deML"
    # PHENIQS_CMD="pheniqs"
    # SAMTOOLS_CMD="samtools"
    # SIMULATE_PY_CMD="$PHENIQS_CODE_HOME/tool/simulate.py"
    # SUMMARIZE_PY_CMD="$PHENIQS_CODE_HOME/tool/summarize.py"
    # EXTRACT_PRIOR_PY_CMD="$PHENIQS_CODE_HOME/tool/prior.py"
    # ANALYZE_PY_CMD="$PHENIQS_CODE_HOME/tool/analyze.py"
    # TODEML_PY_CMD="$PHENIQS_CODE_HOME/tool/todeml.py"
    # SIMULATION_CONFIGURATION="$PHENIQS_CODE_HOME/simulation/core.json"
    #
    # PHENIQS_DEMUX_CONFIG="$PHENIQS_CODE_HOME/simulation/$SIMULATION_PRESET/pheniqs_annotate.json"
    # DEML_INDEX_FILE="$PHENIQS_CODE_HOME/simulation/$SIMULATION_PRESET/deml_index.txt"
    # SOURCE_DATA="~/Downloads/A5KVK_tiny.bam"
    #
    #
    # make_simulated_for_pheniqs() {
    #     echo "\
    #     time $SAMTOOLS_CMD view $SOURCE_DATA | $SIMULATE_PY_CMD \
    #     --configuration $SIMULATION_CONFIGURATION \
    #     --preset $1 \
    #     --report $1_model.json | \
    #     $SAMTOOLS_CMD view -b > $1_simulated.bam \
    #     "
    # }
    #
    # make_simulated_for_deml() {
    #     echo "\
    #     time $SAMTOOLS_CMD view $1_simulated.bam | $TODEML_PY_CMD | \
    #     $SAMTOOLS_CMD view -b > $1_simulated_deml_syntax.bam \
    #     "
    # }
    #
    # pheniqs_sense_prior() {
    #     echo "\
    #     time $PHENIQS_CMD mux \
    #     --sense-input \
    #     --output /dev/null \
    #     --config $PHENIQS_DEMUX_CONFIG \
    #     --input $1_simulated.bam \
    #     --report $1_uniform_report.json \
    #     "
    # }
    #
    # pheniqs_adjust_prior() {
    #     echo "\
    #     time $EXTRACT_PRIOR_PY_CMD \
    #     --original $PHENIQS_DEMUX_CONFIG \
    #     --model $1_model.json \
    #     --report $1_uniform_report.json \
    #     --adjusted $1_adjusted.json \
    #     "
    # }
    #
    # pheniqs_demux() {
    #     echo "\
    #     time $PHENIQS_CMD mux \
    #     --sense-input \
    #     --input $1_simulated.bam \
    #     --config $1_adjusted.json \
    #     --output $1_pheniqs_demux.bam \
    #     --report $1_pheniqs_demux_report.json \
    #     "
    # }
    #
    # deml_demux() {
    #     echo "\
    #     time $DEML_CMD \
    #     --index $DEML_INDEX_FILE \
    #     --outfile $1_deml_demux.bam \
    #     --summary $1_deml_summary.txt \
    #     $1_simulated_deml_syntax.bam \
    #     "
    # }
    #
    # pheniqs_analyze() {
    #     echo "\
    #     time $SAMTOOLS_CMD view $1_pheniqs_demux.bam | \
    #     $ANALYZE_PY_CMD --report $1_pheniqs_analysis.json \
    #     $1_model.json \
    #     "
    # }
    #
    # deml_analyze() {
    #     echo "\
    #     time $SAMTOOLS_CMD view $1_deml_demux.bam | \
    #     $ANALYZE_PY_CMD --report $1_deml_analysis.json \
    #     $1_model.json \
    #     "
    # }
    #
    # summarize() {
    #     echo "\
    #     time $SUMMARIZE_PY_CMD . \
    #     "
    # }
    #
    # do_good() {
    #     make_simulated_for_pheniqs $1;
    #     make_simulated_for_deml $1;
    #     pheniqs_sense_prior $1;
    #     pheniqs_adjust_prior $1;
    #     pheniqs_demux $1;
    #     pheniqs_analyze $1;
    #     deml_demux $1;
    #     deml_analyze $1;
    #     summarize $1;
    # }
    #
    # do_good $SIMULATION_PRESET

from core import *
from simulate_barcode import SimulateBarcode
from simulate_error import SimulateError
from todeml import ToDeML
from prior import Prior
from analyze import Analyze
from collect import Collect
from summarize import Summarize

class Benchmark(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)

    def execute(self):
        if self.action == 'simulate_barcode':
            self.simulate_barcode(self.ontology)

        elif self.action == 'simulate_error':
            self.simulate_error(self.ontology)

        elif self.action == 'todeml':
            self.todeml(self.ontology)

        elif self.action == 'sense_prior':
            self.sense_prior(self.ontology)

        elif self.action == 'adjust_prior':
            self.adjust_prior(self.ontology)

        elif self.action == 'demux_pheniqs':
            self.demux_pheniqs(self.ontology)

        elif self.action == 'demux_deml':
            self.demux_deml(self.ontology)

        elif self.action == 'analyze_pheniqs':
            self.analyze_pheniqs(self.ontology)

        elif self.action == 'analyze_deml':
            self.analyze_deml(self.ontology)

        elif self.action == 'collect':
            self.collect(self.ontology)

        elif self.action == 'summarize':
            self.summarize(self.ontology)

    def simulate_barcode(self, ontology):
        self.log.info('simulating barcode indices')
        pipeline = SimulateBarcode(ontology)
        pipeline.execute()

    def simulate_error(self, ontology):
        self.log.info('simulating errors on barcode indices')
        pipeline = SimulateError(ontology)
        pipeline.execute()

    def todeml(self, ontology):
        self.log.info('transcoding simulated data to deML syntax')
        pipeline = ToDeML(ontology)
        pipeline.execute()

    def sense_prior(self, ontology):
        self.log.info('estimating priors')
        command = [ 'pheniqs', 'mux', '--sense-input', '--output', '/dev/null' ]

        self.log.info('estimating priors for %s with configuration %s', self.instruction['input'], self.instruction['configuration'])

        command.append('--config')
        command.append(self.instruction['configuration'])

        command.append('--input')
        command.append(self.instruction['input'])

        command.append('--report')
        command.append(self.instruction['report'])

        self.log.debug(' '.join([str(i) for i in command]))
        process = Popen(
            args=command,
            # env=self.env,
            cwd=self.current_working_directoy,
            # stdout=PIPE,
            # stderr=PIPE
        )
        output, error = process.communicate()
        code = process.returncode
        if code == 0:
            # self.stdout.write(output.decode('utf8'))
            # self.stderr.write(error.decode('utf8'))
            pass
        else:
            print(output.decode('utf8'))
            print(error.decode('utf8'))
            raise CommandFailedError('pheniqs returned {} when estimating prior'.format(code))

    def adjust_prior(self, ontology):
        self.log.info('adjusting priors')
        pipeline = Prior(ontology)
        pipeline.execute()

    def demux_pheniqs(self, ontology):
        self.log.info('demultiplexing with pheniqs')
        command = [ 'pheniqs', 'mux', '--sense-input' ]

        self.log.info('pheniqs demultiplexing %s with configuration %s', self.instruction['input'], self.instruction['configuration'])

        command.append('--config')
        command.append(self.instruction['configuration'])

        command.append('--input')
        command.append(self.instruction['input'])

        command.append('--output')
        command.append(self.instruction['output'])

        command.append('--report')
        command.append(self.instruction['report'])

        self.log.debug(' '.join([str(i) for i in command]))
        process = Popen(
            args=command,
            # env=self.env,
            cwd=self.current_working_directoy,
            # stdout=PIPE,
            # stderr=PIPE
        )
        output, error = process.communicate()
        code = process.returncode
        if code == 0:
            # self.stdout.write(output.decode('utf8'))
            # self.stderr.write(error.decode('utf8'))
            pass
        else:
            print(output.decode('utf8'))
            print(error.decode('utf8'))
            raise CommandFailedError('pheniqs returned {} when demultiplexing'.format(code))

    def demux_deml(self, ontology):
        self.log.info('demultiplexing with deML')
        command = [ 'deML' ]

        self.log.info('deML demultiplexing %s with configuration %s', self.instruction['input'], self.instruction['configuration'])

        command.append('--index')
        command.append(self.instruction['configuration'])

        command.append('--outfile')
        command.append(self.instruction['output'])

        command.append('--summary')
        command.append(self.instruction['report'])

        command.append(self.instruction['input'])

        self.log.debug(' '.join([str(i) for i in command]))
        process = Popen(
            args=command,
            # env=self.env,
            cwd=self.current_working_directoy,
            # stdout=PIPE,
            # stderr=PIPE
        )
        output, error = process.communicate()
        code = process.returncode
        if code == 0:
            # self.stdout.write(output.decode('utf8'))
            # self.stderr.write(error.decode('utf8'))
            pass
        else:
            print(output.decode('utf8'))
            print(error.decode('utf8'))
            raise CommandFailedError('deML returned {} when demultiplexing'.format(code))

    def analyze_pheniqs(self, ontology):
        self.log.info('analyzing pheniqs results')
        ontology['instruction']['tool'] = 'pheniqs'
        pipeline = Analyze(ontology)
        pipeline.execute()

    def analyze_deml(self, ontology):
        self.log.info('analyzing deML results')
        ontology['instruction']['tool'] = 'deml'
        pipeline = Analyze(ontology)
        pipeline.execute()

    def collect(self, ontology):
        pipeline = Collect(ontology)
        pipeline.execute()

    def summarize(self, ontology):
        pipeline = Summarize(ontology)
        pipeline.execute()

def main():
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    pipeline = None

    try:
        command = CommandLineParser('benchmark')
        if command.help_triggered:
            command.help()
            sys.exit(0)
        else:
            if 'verbosity' in command.instruction and command.instruction['verbosity']:
                logging.getLogger().setLevel(log_levels[command.instruction['verbosity']])

            pipeline = Benchmark(command.configuration)
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
