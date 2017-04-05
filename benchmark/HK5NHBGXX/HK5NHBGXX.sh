#!/bin/zsh

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

LG_HOME="/home/lg1883"
PROCESSORS=16
INPUT_BASE="${LG_HOME}/flowcell"
PHENIQS_HOME="${LG_HOME}/code/pheniqs"
FASQ_MULTX="${LG_HOME}/code/fastq-multx/fastq-multx"
PHENIQS="$PHENIQS_HOME/pheniqs"
FLOWCELL_ID="HK5NHBGXX"
CONFIG_FOLDER="$PHENIQS_HOME/benchmark/${FLOWCELL_ID}"
FASTQ_MULTX_BARCODES="${CONFIG_FOLDER}/multx.tsv"
BENCHMARK_FOLDER="${LG_HOME}/benchmark/${FLOWCELL_ID}"
LOG_FILE="$BENCHMARK_FOLDER/benchmark.log"

function clear_os_cache() {
    sync;
    echo 1 > /proc/sys/vm/drop_caches;
    echo 2 > /proc/sys/vm/drop_caches;
    echo 3 > /proc/sys/vm/drop_caches;
    echo 1 > /proc/sys/vm/compact_memory;
};

function run_fastq_multx_split() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/multx/split";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark fastq-multx mdd split fastq" >> $LOG_FILE

    clear_os_cache
    {   time $FASQ_MULTX -B "$FASTQ_MULTX_BARCODES" \
        -d 1 \
        -q 1 \
        ${INPUT_BASE}/HK5NHBGXX_l01n02.33200000038101.fastq.gz \
        ${INPUT_BASE}/HK5NHBGXX_l01n03.3330000003810e.fastq.gz \
        ${INPUT_BASE}/HK5NHBGXX_l01n01.33100000038104.fastq.gz \
        ${INPUT_BASE}/HK5NHBGXX_l01n04.3340000003810b.fastq.gz \
        -o n/a \
        -o n/a \
        -o "$OUTPUT_FOLDER/%_1.fastq.gz" \
        -o "$OUTPUT_FOLDER/%_2.fastq.gz" > /dev/null 2>&1; \
    } 2>> $LOG_FILE
};

function run_pheniqs_mdd_fastq_split() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/mdd/fastq/split";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs mdd split fastq" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/fastq_split.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder mdd \
        --quality \
        > /dev/null; \
    } 2>> $LOG_FILE
};
function run_pheniqs_mdd_fastq_split_quality() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/mdd/fastq/split";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs mdd split fastq quality" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/fastq_split.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder mdd \
        > $OUTPUT_FOLDER/report.json; \
    } 2>> $LOG_FILE
};

function run_pheniqs_mdd_fastq_interleaved() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/mdd/fastq/interleaved";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs mdd interleaved fastq" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/fastq_interleaved.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder mdd \
        --quality \
        > /dev/null; \
    } 2>> $LOG_FILE
};
function run_pheniqs_mdd_fastq_interleaved_quality() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/mdd/fastq/interleaved";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs mdd interleaved fastq quality" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/fastq_interleaved.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder mdd \
        > $OUTPUT_FOLDER/report.json; \
    } 2>> $LOG_FILE
};

function run_pheniqs_mdd_cram_interleaved() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/mdd/cram/interleaved";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs mdd interleaved cram" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/cram_interleaved.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder mdd \
        --quality \
        > /dev/null; \
    } 2>> $LOG_FILE
};
function run_pheniqs_mdd_cram_interleaved_quality() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/mdd/cram/interleaved";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs mdd interleaved cram quality" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/cram_interleaved.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder mdd \
        > $OUTPUT_FOLDER/report.json; \
    } 2>> $LOG_FILE
};

function run_pheniqs_mdd_cram_combined() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/mdd/cram/combined";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs mdd combined cram" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/cram_combined.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder mdd \
        --quality \
        > /dev/null; \
    } 2>> $LOG_FILE
};
function run_pheniqs_mdd_cram_combined_quality() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/mdd/cram/combined";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs mdd combined cram quality" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/cram_combined.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder mdd \
        > $OUTPUT_FOLDER/report.json; \
    } 2>> $LOG_FILE
};

function run_pheniqs_pamld_fastq_split() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/pamld/cram/interleaved";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs pamld split fastq" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/fastq_split.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder pamld \
        --quality \
        > /dev/null; \
    } 2>> $LOG_FILE
};
function run_pheniqs_pamld_fastq_split_quality() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/pamld/cram/interleaved";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs pamld split fastq quality" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/fastq_split.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder pamld \
        > $OUTPUT_FOLDER/report.json; \
    } 2>> $LOG_FILE
};

function run_pheniqs_pamld_fastq_interleaved() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/pamld/fastq/interleaved";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs pamld interleaved fastq" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/fastq_interleaved.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder pamld \
        --quality \
        > /dev/null; \
    } 2>> $LOG_FILE
};
function run_pheniqs_pamld_fastq_interleaved_quality() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/pamld/fastq/interleaved";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs pamld interleaved fastq quality" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/fastq_interleaved.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder pamld \
        > $OUTPUT_FOLDER/report.json; \
    } 2>> $LOG_FILE
};

function run_pheniqs_pamld_cram_interleaved() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/pamld/fastq/split";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs pamld interleaved cram" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/cram_interleaved.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder pamld \
        --quality \
        > /dev/null; \
    } 2>> $LOG_FILE
};
function run_pheniqs_pamld_cram_interleaved_quality() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/pamld/fastq/split";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs pamld interleaved cram quality" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/cram_interleaved.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder pamld \
        > $OUTPUT_FOLDER/report.json; \
    } 2>> $LOG_FILE
};

function run_pheniqs_pamld_cram_combined() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/pamld/cram/combined";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs pamld combined cram" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/cram_combined.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder pamld \
        --quality \
        > /dev/null; \
    } 2>> $LOG_FILE
};
function run_pheniqs_pamld_cram_combined_quality() {
    OUTPUT_FOLDER="$BENCHMARK_FOLDER/pheniqs/pamld/cram/combined";
    mkdir -p "$OUTPUT_FOLDER";

    echo -e "\n${FLOWCELL_ID} benchmark pheniqs pamld combined cram quality" >> $LOG_FILE

    clear_os_cache
    {   time $PHENIQS demux \
        --config "${CONFIG_FOLDER}/cram_combined.json" \
        --base-output "$OUTPUT_FOLDER" \
        --threads ${PROCESSORS} \
        --decoder pamld \
        > $OUTPUT_FOLDER/report.json; \
    } 2>> $LOG_FILE
};

rm -rf "$BENCHMARK_FOLDER"
mkdir -p "$BENCHMARK_FOLDER"

run_pheniqs_mdd_fastq_split
run_pheniqs_mdd_fastq_split_quality

run_pheniqs_pamld_fastq_split
run_pheniqs_pamld_fastq_split_quality

run_pheniqs_mdd_fastq_interleaved
run_pheniqs_mdd_fastq_interleaved_quality

run_pheniqs_pamld_fastq_interleaved
run_pheniqs_pamld_fastq_interleaved_quality

run_pheniqs_mdd_cram_interleaved
run_pheniqs_mdd_cram_interleaved_quality

run_pheniqs_pamld_cram_interleaved
run_pheniqs_pamld_cram_interleaved_quality

run_pheniqs_mdd_cram_combined
run_pheniqs_mdd_cram_combined_quality

run_pheniqs_pamld_cram_combined
run_pheniqs_pamld_cram_combined_quality

run_fastq_multx_split
