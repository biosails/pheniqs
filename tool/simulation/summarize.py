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

import logging
from core.error import *
from core import Job
from core import to_json

def fscore(precision, recall):
    if(precision + recall) > 0:
        return 2.0 * (precision * recall) / (precision + recall)
    else:
        return 0

class Summarize(Job):
    def __init__(self, ontology):
        Job.__init__(self, ontology)
        self.log = logging.getLogger('Summarize')

        if self.instruction['refresh'] and 'summary' in self.barcode_simulation:
            del self.barcode_simulation['summary']

    @property
    def bsid(self):
        return self.barcode_simulation['barcode']['model']['genealogy']['bsid']

    @property
    def barcode_simulation(self):
        return self.ontology['barcode simulation']

    @property
    def barcode_simulation_summary(self):
        if 'summary' not in self.session:
            self.session['summary'] = self.summarize_simulation()
        return self.session['summary']

    @property
    def barcode_binning_model(self):
        if 'barcode binning model' not in self.session:
            self.session['barcode binning model'] = self.compose_barcode_binning_model()
        return self.session['barcode binning model']

    def initialize_accurecy_record(self, record=None):
        for key, value in {
            'TP': 0,
            'FP': 0,
            'FN': 0,
            'TP_FN': 0,
            'TP_FP': 0,
            'FDR': 0,
            'MR': 0,
            'precision': 0,
            'recall': 0,
            'fscore': 0,
        }.items(): record[key] = value

    def finalize_accurecy_record(self, record):
        record['TP_FP'] = record['TP'] + record['FP']
        if record['TP_FP'] > 0:
            record['precision']    = record['TP'] / record['TP_FP']
            record['FDR']          = record['FP'] / record['TP_FP']

        record['TP_FN'] = record['TP'] + record['FN']
        if record['TP_FN'] > 0:
            record['recall']       = record['TP'] / record['TP_FN']
            record['MR']           = record['FN'] / record['TP_FN']

        record['fscore'] = fscore(record['precision'], record['recall'])

    def compose_barcode_binning_model(self):
        def assemble_classifier_binning_map(classifier_model, classifier_binning_model):
            if 'codec' in classifier_model:
                classifier_binning_model['barcode bin by index'] = []
                for bin_index, bin_limit in enumerate(classifier_binning_model['interval']):
                    classifier_binning_model['barcode bin by index'].append({
                        'index': bin_index,
                        'limit': bin_limit,
                        'barcode': [],
                    })

                barcode_density_stack = []
                for barcode in classifier_model['codec'].values():
                    barcode_density_stack.append(
                        {
                            'index': barcode['index'],
                            'count': barcode['count'],
                            'density': barcode['count'] / classifier_binning_model['count'],
                        }
                    )
                barcode_density_stack.sort(key=lambda i: i['density'])

                order = 0
                while barcode_density_stack:
                    record = barcode_density_stack.pop(0)
                    record['order'] = order
                    order += 1
                    for bin_model in classifier_binning_model['barcode bin by index']:
                        if record['density'] < bin_model['limit']:
                            record['bin index'] = bin_model['index']
                            bin_model['barcode'].append(record)
                            break

        barcode_model = self.barcode_simulation['barcode']['model']
        barcode_binning_model = {}
        interval = [
            0.001,
            0.003,
            0.01,
            0.03,
            1.0,
        ]
        for classifier_type in [ 'multiplex', 'cellular', 'molecular' ]:
            if classifier_type in barcode_model:
                classifier_model = barcode_model[classifier_type]
                if isinstance(classifier_model, dict):
                    classifier_binning_model = {
                        'index': classifier_model['index'],
                        'count': barcode_model['count'],
                        'classifier type': classifier_type,
                        'interval': interval,
                    }
                    assemble_classifier_binning_map(classifier_model, classifier_binning_model)
                    barcode_binning_model[classifier_type] = classifier_binning_model

                elif isinstance(classifier_model, list):
                    barcode_binning_model[classifier_type] = []
                    for item in classifier_model:
                        classifier_binning_model = {
                            'index': item['index'],
                            'count': barcode_model['count'],
                            'classifier type': classifier_type,
                            'interval': interval,
                        }
                        assemble_classifier_binning_map(item, classifier_binning_model)
                        barcode_binning_model[classifier_type].append(classifier_binning_model)

        return barcode_binning_model

    def summarize_simulation(self):
        def summarize_classifier_simulation(classifier_analysis, classifier_report):
            # classification code from analysis
            # if true_barcode['index'] > 0:
            #     if true_barcode['index'] == decoded_barcode['index']:
            #         true_barcode['accumulate']['real'][qc]['TP'] += 1
            #     else:
            #         if decoded_barcode['index'] > 0:
            #             true_barcode['accumulate']['real'][qc]['FN'] += 1
            #             decoded_barcode['accumulate']['real'][qc]['FP'] += 1
            #         else:
            #             true_barcode['accumulate']['noise'][qc]['FN'] += 1
            #             decoded_barcode['accumulate']['noise'][qc]['FP'] += 1 This is actually a real read
            # else:
            #     if decoded_barcode['index'] > 0:
            #         true_barcode['accumulate']['noise'][qc]['FN'] += 1
            #         decoded_barcode['accumulate']['noise'][qc]['FP'] += 1
            #     else:
            #         true_barcode['accumulate']['noise'][qc]['TP'] += 1
            classifier_report['classified count'] = classifier_analysis['classified count']
            classifier_report['rate'] = classifier_analysis['simulated substitution rate']
            classifier_report['simulated rate'] = classifier_analysis['simulated substitution rate']
            classifier_report['expected rate'] = classifier_analysis['expected substitution rate']
            if 'requested substitution rate' in classifier_analysis:
                classifier_report['requested rate'] = classifier_analysis['requested substitution rate']
            else:
                classifier_report['requested rate'] = 0

            classifier_report['barcode'] = []
            classifier_report['classified'] = {
                'count': 0,
            }
            self.initialize_accurecy_record(classifier_report['classified'])

            classifier_report['unclassified'] = {
                'count': None,
                'index': None,
            }
            self.initialize_accurecy_record(classifier_report['unclassified'])

            classifier_report['summary'] = {
                'count': 0,
            }
            self.initialize_accurecy_record(classifier_report['summary'])

            if 'codec' in classifier_analysis:
                codec = classifier_analysis['codec']
                unclassified_report = classifier_report['unclassified']
                for barcode_analysis in codec.values():
                    if 'accumulate' in barcode_analysis:
                        # real,pass,TP  :   Real read was correctly classified and passing filter.
                        # real,fail,TP  :   Real read was correctly classified but erroneously filtered.
                        # real,pass,FP  :   Cross contamination erroneously passing filter.
                        # real,fail,FP  :   Cross contamination filtered.
                        # noise,pass,FP :   Noise erroneously decoded as real and erroneously passing filter
                        # noise,fail,FP :   Noise erroneously decoded as real and filtered
                        # real,pass,FN  :   Cross contamination erroneously passing filter.
                        # real,fail,FN  :   Cross contamination filtered.
                        # noise,pass,FN :   Real erroneously decoded as noise. Only in MDD marked as pass.
                        # noise,fail,FN :   Real erroneously decoded as noise and filtered.
                        accumulate = barcode_analysis['accumulate']
                        barcode_report = {
                            'index': barcode_analysis['index'],
                            'count': barcode_analysis['count'],
                        }
                        self.initialize_accurecy_record(barcode_report)
                        classifier_report['barcode'].append(barcode_report)

                        # Real read was correctly classified and passing filter.
                        barcode_report['TP'] += accumulate['real']['pass']['TP']

                        # Real read was correctly classified but erroneously filtered.
                        barcode_report['FN'] += accumulate['real']['fail']['TP']
                        unclassified_report['FP'] += accumulate['real']['fail']['TP']

                        # Cross contamination erroneously passing filter.
                        barcode_report['FP'] += accumulate['real']['pass']['FP']

                        # Cross contamination filtered.
                        unclassified_report['FP'] += accumulate['real']['fail']['FP']

                        # Noise erroneously decoded as real and erroneously passing filter
                        barcode_report['FP'] += accumulate['noise']['pass']['FP']

                        # Noise erroneously decoded as real and filtered
                        # already accounted for with ['noise']['fail']['FN'] on the noise bin
                        # unclassified_analysis['TP'] += accumulate['noise']['fail']['FP']

                        # Cross contamination erroneously passing filter
                        barcode_report['FN'] += accumulate['real']['pass']['FN']

                        # Cross contamination filtered.
                        barcode_report['FN'] += accumulate['real']['fail']['FN']

                        # Real erroneously decoded as noise. Only in MDD marked as pass.
                        barcode_report['FN'] += accumulate['noise']['pass']['FN']

                        # Real erroneously decoded as noise and filtered.
                        barcode_report['FN'] += accumulate['noise']['fail']['FN']

                        classifier_report['classified']['count'] += barcode_report['count']
                        classifier_report['classified']['TP'] += barcode_report['TP']
                        classifier_report['classified']['FP'] += barcode_report['FP']
                        classifier_report['classified']['FN'] += barcode_report['FN']
                        self.finalize_accurecy_record(barcode_report)
                        self.log.debug (
                            '%s %s %-9s Index %-5s %-9s %-9s %s',
                            ssid,
                            classifier_type,
                            tool,
                            barcode_report['index'],
                            barcode_report['count'],
                            barcode_report['TP_FN'],
                            barcode_report['count'] - barcode_report['TP_FN'],
                        )

                classifier_report['summary']['count'] += classifier_report['classified']['count']
                classifier_report['summary']['TP'] += classifier_report['classified']['TP']
                classifier_report['summary']['FP'] += classifier_report['classified']['FP']
                classifier_report['summary']['FN'] += classifier_report['classified']['FN']

            if 'unclassified' in classifier_analysis:
                barcode_analysis = classifier_analysis['unclassified']
                if 'accumulate' in barcode_analysis:
                    # noise,pass,TP :   Noise read correctly classified. Only in MDD marked as pass.
                    # noise,fail,TP :   Noise read correctly classified and filtered.
                    # noise,pass,FP :   Real erroneously decoded as noise. Only in MDD marked as pass.
                    # noise,fail,FP :   Real erroneously decoded as noise and filtered.
                    # noise,pass,FN :   Noise erroneously decoded as real and erroneously passing filter
                    # noise,fail,FN :   Noise erroneously decoded as real and filtered
                    accumulate = barcode_analysis['accumulate']
                    barcode_report = classifier_report['unclassified']
                    barcode_report['index'] = barcode_analysis['index']
                    barcode_report['count'] = barcode_analysis['count']

                    # Noise read correctly classified. Only in MDD marked as pass.
                    barcode_report['TP'] += accumulate['noise']['pass']['TP']

                    # Noise read correctly classified and filtered.
                    barcode_report['TP'] += accumulate['noise']['fail']['TP']

                    # Real erroneously decoded as noise. Only in MDD marked as pass.
                    barcode_report['FP'] += accumulate['noise']['pass']['FP']

                    # Real erroneously decoded as noise and filtered.
                    barcode_report['FP'] += accumulate['noise']['fail']['FP']

                    # Noise erroneously decoded as real and erroneously passing filter
                    barcode_report['FN'] += accumulate['noise']['pass']['FN']

                    # Noise erroneously decoded as real and filtered
                    barcode_report['TP'] += accumulate['noise']['fail']['FN']

                    classifier_report['summary']['count'] += barcode_report['count']
                    classifier_report['summary']['TP'] += barcode_report['TP']
                    classifier_report['summary']['FP'] += barcode_report['FP']
                    classifier_report['summary']['FN'] += barcode_report['FN']
                    self.finalize_accurecy_record(barcode_report)
                    self.log.debug (
                        '%s %s %-9s %-10s %-9s %-9s %s',
                        ssid,
                        classifier_type,
                        tool,
                        'Noise ',
                        barcode_report['count'],
                        barcode_report['TP_FN'],
                        barcode_report['count'] - barcode_report['TP_FN'],
                    )

            self.finalize_accurecy_record(classifier_report['classified'])
            self.finalize_accurecy_record(classifier_report['summary'])
            self.log.debug (
                '%s %s %-9s Balance    %-9s %-9s %s',
                ssid,
                classifier_type,
                tool,
                classifier_report['count'],
                classifier_report['summary']['TP_FN'],
                classifier_report['count'] - classifier_report['summary']['TP_FN'],
            )
            self.log.debug (
                '%s:%s:%-9s Classified  %-9s %-9s %s',
                ssid,
                classifier_type,
                tool,
                classifier_report['classified count'],
                classifier_report['classified']['TP_FN'],
                classifier_report['classified count'] - classifier_report['classified']['TP_FN'],
            )

        self.log.info('summarizing barcode simulation %s', self.bsid)
        barcode_simulation_summary = {}
        for ssid, substitution_analysis in self.barcode_simulation['substitution'].items():
            substitution_report = {
                'ssid': ssid,
                'count': substitution_analysis['model']['count'],
                'tool': {}
            }
            for tool, tool_analysis_key in {
                'deml': 'deml demultiplex analysis',
                'pamld': 'pamld demultiplex analysis',
                'pamld_ap': 'pamld accurate prior demultiplex analysis',
                'pamld_u': 'pamld uniform demultiplex analysis',
                'mdd': 'mdd demultiplex analysis',
            }.items():
                if tool_analysis_key in substitution_analysis:
                    tool_substitution_analysis = substitution_analysis[tool_analysis_key]
                    tool_report = {
                        'tool': tool,
                        'classifier': {}
                    }
                    substitution_report['tool'][tool] = tool_report

                    for classifier_type in [ 'multiplex', 'cellular', 'molecular' ]:
                        if classifier_type in tool_substitution_analysis:
                            classifier_analysis = tool_substitution_analysis[classifier_type]
                            if isinstance(classifier_analysis, dict):
                                classifier_report = {
                                    'index': classifier_analysis['index'],
                                    'classifier type': classifier_type,
                                    'count': substitution_report['count'],
                                }
                                summarize_classifier_simulation(classifier_analysis, classifier_report)
                                tool_report['classifier'][classifier_type] = classifier_report

                            elif isinstance(classifier_analysis, list):
                                tool_report['classifier'][classifier_type] = []
                                for item in classifier_analysis:
                                    classifier_report = {
                                        'index': item['index'],
                                        'classifier type': classifier_type,
                                        'count': substitution_report['count'],
                                    }
                                    summarize_classifier_simulation(item, classifier_report)
                                    tool_report['classifier'][classifier_type].append(classifier_report)

            if ssid not in barcode_simulation_summary:
                barcode_simulation_summary[ssid] = substitution_report
            else:
                self.log.warning('substitution %s already present', ssid)

        self.summarize_binned_simulation(barcode_simulation_summary)
        return barcode_simulation_summary

    def summarize_binned_simulation(self, barcode_simulation_summary):
        def assemble_classifier_binned_barcode_simulation_summary(classifier_binning_model, classifier_report):
            binned_classifier_report = []
            barcode_report_by_index = {}
            for barcode_report in classifier_report['barcode']:
                barcode_report_by_index[barcode_report['index']] = barcode_report

            for bin_model in classifier_binning_model['barcode bin by index']:
                barcode_bin_report = {
                    'index': bin_model['index'],
                    'limit': bin_model['limit'],
                    'barcode count': 0,
                }
                self.initialize_accurecy_record(barcode_bin_report)
                for barcode in bin_model['barcode']:
                    barcode_report = barcode_report_by_index[barcode['index']]
                    barcode_bin_report['TP'] += barcode_report['TP']
                    barcode_bin_report['FP'] += barcode_report['FP']
                    barcode_bin_report['FN'] += barcode_report['FN']
                    barcode_bin_report['barcode count'] += 1
                    self.finalize_accurecy_record(barcode_bin_report)
                binned_classifier_report.append(barcode_bin_report)

            classifier_report['binned barcode'] = binned_classifier_report

        self.log.info('summarizing binned barcode simulation %s', self.bsid)
        for ssid, substitution_report in barcode_simulation_summary.items():
            for tool, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    if classifier_type in self.barcode_binning_model:
                        classifier_binning_model = self.barcode_binning_model[classifier_type]

                        if isinstance(classifier_report, dict):
                            assemble_classifier_binned_barcode_simulation_summary(classifier_binning_model, classifier_report)

                        elif isinstance(classifier_report, list):
                            for model, report in zip(classifier_binning_model, classifier_report):
                                assemble_classifier_binned_barcode_simulation_summary(model, report)

    def execute(self):
        self.log.info('summarizing %s benchmarks', self.instruction['preset'])

        if self.instruction['preset'] == 'json':
            print(to_json(self.barcode_simulation_summary))

        elif self.instruction['preset'] == 'bin_model':
            print(to_json(self.barcode_binning_model))

        elif self.instruction['preset'] == 'decoder_summary_R':
            self.summarize_decoder_accuracy_benchmark_R()

        elif self.instruction['preset'] == 'decoder_summary':
            self.summarize_decoder_accuracy_benchmark()

        elif self.instruction['preset'] == 'barcode_summary_R':
            self.summarize_barcode_accuracy_benchmark_R()

        elif self.instruction['preset'] == 'barcode_summary':
            self.summarize_barcode_accuracy_benchmark()

        elif self.instruction['preset'] == 'noise_summary_R':
            self.summarize_noise_accuracy_benchmark_R()

        elif self.instruction['preset'] == 'noise_summary':
            self.summarize_noise_accuracy_benchmark()

        elif self.instruction['preset'] == 'classified_summary_R':
            self.summarize_classified_accuracy_benchmark_R()

        elif self.instruction['preset'] == 'classified_summary':
            self.summarize_classified_accuracy_benchmark()

        elif self.instruction['preset'] == 'binned_decoder_summary_R':
            self.summarize_binned_decoder_accuracy_benchmark_R()

        elif self.instruction['preset'] == 'binned_decoder_summary':
            self.summarize_binned_decoder_accuracy_benchmark()

        elif self.instruction['preset'] == 'quality_distribution':
            self.summarize_decoder_quality_distribution_R()

        elif self.instruction['preset'] == 'barcode_distribution':
            self.summarize_barcode_distribution_R()

        self.close()

    def summarize_decoder_accuracy_benchmark_R(self):
        header = [
            'ssid',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'variable',
            'value',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    for variable in [ 'FDR', 'MR', 'TP', 'FN', 'FP', 'TP_FP', 'TP_FN', 'precision', 'recall', 'fscore' ]:
                        row = [
                            ssid,
                            classifier_report['simulated rate'],
                            classifier_report['expected rate'],
                            classifier_report['requested rate'],
                            tool_id,
                            classifier_type,
                            variable,
                            classifier_report['summary'][variable],
                        ]
                        table.append(row)

        table.sort(key=lambda i: i[6])
        table.sort(key=lambda i: i[5])
        table.sort(key=lambda i: i[4])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_decoder_accuracy_benchmark(self):
        header = [
            'ssid',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'TP',
            'FP',
            'FN',
            'TP_FN',
            'TP_FP',
            'FDR',
            'MR',
            'precision',
            'recall',
            'fscore',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    row = [
                        ssid,
                        classifier_report['simulated rate'],
                        classifier_report['expected rate'],
                        classifier_report['requested rate'],
                        tool_id,
                        classifier_type,
                        classifier_report['summary']['TP'],
                        classifier_report['summary']['FP'],
                        classifier_report['summary']['FN'],
                        classifier_report['summary']['TP_FN'],
                        classifier_report['summary']['TP_FP'],
                        classifier_report['summary']['FDR'],
                        classifier_report['summary']['MR'],
                        classifier_report['summary']['precision'],
                        classifier_report['summary']['recall'],
                        classifier_report['summary']['fscore'],
                    ]
                    table.append(row)

        table.sort(key=lambda i: i[5])
        table.sort(key=lambda i: i[4])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_barcode_accuracy_benchmark_R(self):
        header = [
            'ssid',
            'index',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'variable',
            'value',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    for barcode_report in classifier_report['barcode']:
                        for variable in [ 'FDR', 'MR', 'TP', 'FN', 'FP', 'TP_FP', 'TP_FN', 'precision', 'recall', 'fscore' ]:
                            row = [
                                ssid,
                                barcode_report['index'],
                                classifier_report['simulated rate'],
                                classifier_report['expected rate'],
                                classifier_report['requested rate'],
                                tool_id,
                                classifier_type,
                                variable,
                                barcode_report[variable],
                            ]
                            table.append(row)

        table.sort(key=lambda i: i[7])
        table.sort(key=lambda i: i[6])
        table.sort(key=lambda i: i[5])
        table.sort(key=lambda i: i[4])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_barcode_accuracy_benchmark(self):
        header = [
            'ssid',
            'index',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'TP',
            'FP',
            'FN',
            'TP_FN',
            'TP_FP',
            'FDR',
            'MR',
            'precision',
            'recall',
            'fscore',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    for barcode_report in classifier_report['barcode']:
                        row = [
                            ssid,
                            barcode_report['index'],
                            classifier_report['simulated rate'],
                            classifier_report['expected rate'],
                            classifier_report['requested rate'],
                            tool_id,
                            classifier_type,
                            barcode_report['TP'],
                            barcode_report['FP'],
                            barcode_report['FN'],
                            barcode_report['TP_FN'],
                            barcode_report['TP_FP'],
                            barcode_report['FDR'],
                            barcode_report['MR'],
                            barcode_report['precision'],
                            barcode_report['recall'],
                            barcode_report['fscore'],
                        ]
                        table.append(row)

        table.sort(key=lambda i: i[6])
        table.sort(key=lambda i: i[5])
        table.sort(key=lambda i: i[4])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_noise_accuracy_benchmark_R(self):
        header = [
            'ssid',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'variable',
            'value',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    for variable in [ 'FDR', 'MR', 'TP', 'FN', 'FP', 'TP_FP', 'TP_FN', 'precision', 'recall', 'fscore' ]:
                        row = [
                            ssid,
                            classifier_report['simulated rate'],
                            classifier_report['expected rate'],
                            classifier_report['requested rate'],
                            tool_id,
                            classifier_type,
                            variable,
                            classifier_report['unclassified'][variable],
                        ]
                        table.append(row)

        table.sort(key=lambda i: i[6])
        table.sort(key=lambda i: i[5])
        table.sort(key=lambda i: i[4])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_noise_accuracy_benchmark(self):
        header = [
            'ssid',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'TP',
            'FP',
            'FN',
            'TP_FN',
            'TP_FP',
            'FDR',
            'MR',
            'precision',
            'recall',
            'fscore',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    row = [
                        ssid,
                        classifier_report['simulated rate'],
                        classifier_report['expected rate'],
                        classifier_report['requested rate'],
                        tool_id,
                        classifier_type,
                        classifier_report['unclassified']['TP'],
                        classifier_report['unclassified']['FP'],
                        classifier_report['unclassified']['FN'],
                        classifier_report['unclassified']['TP_FN'],
                        classifier_report['unclassified']['TP_FP'],
                        classifier_report['unclassified']['FDR'],
                        classifier_report['unclassified']['MR'],
                        classifier_report['unclassified']['precision'],
                        classifier_report['unclassified']['recall'],
                        classifier_report['unclassified']['fscore'],
                    ]
                    table.append(row)

        table.sort(key=lambda i: i[5])
        table.sort(key=lambda i: i[4])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_classified_accuracy_benchmark_R(self):
        header = [
            'ssid',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'variable',
            'value',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    for variable in [ 'FDR', 'MR', 'TP', 'FN', 'FP', 'TP_FP', 'TP_FN', 'precision', 'recall', 'fscore' ]:
                        row = [
                            ssid,
                            classifier_report['simulated rate'],
                            classifier_report['expected rate'],
                            classifier_report['requested rate'],
                            tool_id,
                            classifier_type,
                            variable,
                            classifier_report['classified'][variable],
                        ]
                        table.append(row)

        table.sort(key=lambda i: i[6])
        table.sort(key=lambda i: i[5])
        table.sort(key=lambda i: i[4])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_classified_accuracy_benchmark(self):
        header = [
            'ssid',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'TP',
            'FP',
            'FN',
            'TP_FN',
            'TP_FP',
            'FDR',
            'MR',
            'precision',
            'recall',
            'fscore',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    row = [
                        ssid,
                        classifier_report['simulated rate'],
                        classifier_report['expected rate'],
                        classifier_report['requested rate'],
                        tool_id,
                        classifier_type,
                        classifier_report['classified']['TP'],
                        classifier_report['classified']['FP'],
                        classifier_report['classified']['FN'],
                        classifier_report['classified']['TP_FN'],
                        classifier_report['classified']['TP_FP'],
                        classifier_report['classified']['FDR'],
                        classifier_report['classified']['MR'],
                        classifier_report['classified']['precision'],
                        classifier_report['classified']['recall'],
                        classifier_report['classified']['fscore'],
                    ]
                    table.append(row)

        table.sort(key=lambda i: i[5])
        table.sort(key=lambda i: i[4])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_binned_decoder_accuracy_benchmark_R(self):
        header = [
            'ssid',
            'bin',
            'count',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'variable',
            'value',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    for barcode_bin_report in classifier_report['binned barcode']:
                        for variable in [ 'FDR', 'MR', 'TP', 'FN', 'FP', 'TP_FP', 'TP_FN', 'precision', 'recall', 'fscore' ]:
                            row = [
                                ssid,
                                barcode_bin_report['index'],
                                barcode_bin_report['barcode count'],
                                classifier_report['simulated rate'],
                                classifier_report['expected rate'],
                                classifier_report['requested rate'],
                                tool_id,
                                classifier_type,
                                variable,
                                barcode_bin_report[variable],
                            ]
                            table.append(row)

        table.sort(key=lambda i: i[8])
        table.sort(key=lambda i: i[7])
        table.sort(key=lambda i: i[6])
        table.sort(key=lambda i: i[3])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_binned_decoder_accuracy_benchmark(self):
        header = [
            'ssid',
            'bin',
            'count',
            'rate',
            'expected',
            'requested',
            'tool',
            'classifier',
            'TP',
            'FP',
            'FN',
            'TP_FN',
            'TP_FP',
            'FDR',
            'MR',
            'precision',
            'recall',
            'fscore',
        ]
        table = []

        for ssid, substitution_report in self.barcode_simulation_summary.items():
            for tool_id, tool_report in substitution_report['tool'].items():
                for classifier_type, classifier_report in tool_report['classifier'].items():
                    for barcode_bin_report in classifier_report['binned barcode']:
                        row = [
                            ssid,
                            barcode_bin_report['index'],
                            barcode_bin_report['barcode count'],
                            classifier_report['simulated rate'],
                            classifier_report['expected rate'],
                            classifier_report['requested rate'],
                            tool_id,
                            classifier_type,
                            barcode_bin_report['TP'],
                            barcode_bin_report['FP'],
                            barcode_bin_report['FN'],
                            barcode_bin_report['TP_FN'],
                            barcode_bin_report['TP_FP'],
                            barcode_bin_report['FDR'],
                            barcode_bin_report['MR'],
                            barcode_bin_report['precision'],
                            barcode_bin_report['recall'],
                            barcode_bin_report['fscore'],
                        ]
                        table.append(row)

        table.sort(key=lambda i: i[3])
        table.sort(key=lambda i: i[2])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_decoder_quality_distribution_R(self):
        header = [
            'ssid',
            'rate',
            'expected',
            'requested',
            'classifier',
            'quality',
            'density',
        ]
        table = []

        for ssid, substitution_analysis in self.barcode_simulation['substitution'].items():
            if 'model' in substitution_analysis:
                substitution_model = substitution_analysis['model']
                for classifier_type in [ 'multiplex', 'cellular', 'molecular' ]:
                    if classifier_type in substitution_model:
                        classifier_model = substitution_model['multiplex']
                        if 'quality distribution' in classifier_model:
                            simulated_rate = classifier_model['simulated substitution rate']
                            expected_rate = classifier_model['expected substitution rate']
                            if 'requested substitution rate' in classifier_model:
                                requested_rate = classifier_model['requested substitution rate']
                            else:
                                requested_rate = 0

                            for quality, density in enumerate(classifier_model['quality distribution']):
                                row = [
                                    ssid,
                                    simulated_rate,
                                    expected_rate,
                                    requested_rate,
                                    classifier_type,
                                    quality,
                                    density
                                ]
                                table.append(row)

        table.sort(key=lambda i: i[5])
        table.sort(key=lambda i: i[4])
        table.sort(key=lambda i: i[1])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))

    def summarize_barcode_distribution_R(self):
        header = [
            'classifier',
            'index',
            'count',
            'fraction',
            'order',
            'bin',
        ]
        table = []

        for classifier_type, classifier_binning_model in self.barcode_binning_model.items():
            for barcode_binning_model in classifier_binning_model['barcode bin by index']:
                for barcode_model in barcode_binning_model['barcode']:

                        row = [
                            classifier_type,
                            barcode_model['index'],
                            barcode_model['count'],
                            barcode_model['density'],
                            barcode_model['order'],
                            barcode_model['bin index'],
                        ]
                        table.append(row)

        table.sort(key=lambda i: i[2])
        table.sort(key=lambda i: i[1])
        table.sort(key=lambda i: i[0])

        print(','.join(header))
        print('\n'.join([','.join([str(field) for field in row]) for row in table]))
