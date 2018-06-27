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

RUN_DIRECTORY="151023_NB501157_0004_AHK5NHBGXX"
FLOWCELL_DIRECTORY="/volume/babylon/brooklyn/work/flowcell/${RUN_DIRECTORY}"
FLOWCELL_ID="HK5NHBGXX"
BASE_FOLDER="/volume/babylon/brooklyn/work/${RUN_DIRECTORY}"
LOG_FILE="$BASE_FOLDER/benchmark.log"
TMP_FOLDER="${BASE_FOLDER}/tmp"
PICARD_BARCODE_DIR="${BASE_FOLDER}/barcodes"
CONFIG_DIRECTORY="/home/lg1883/code/pheniqs/benchmark/${FLOWCELL_ID}"
PICARD_READ_STRUCTURE="36T8B1S8B1S36T"
function clear_os_cache() {
    sync;
    echo 1 > /proc/sys/vm/drop_caches;
};
function run_picard_ExtractIlluminaBarcodes() {
    OUTPUT_FOLDER="${BASE_FOLDER}/ExtractIlluminaBarcodes"
    BARCODE_FILE_PATH="${CONFIG_DIRECTORY}/picard_barcode_file.tsv"

    rm -rf "$OUTPUT_FOLDER"
    mkdir -p "$OUTPUT_FOLDER"

    clear_os_cache

    echo -e "\nbenchmark ExtractIlluminaBarcodes" >> $LOG_FILE
    {
        time \
        /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx4g \
        -jar /opt/local/tubo/picard-tools-2.7.1/picard.jar \
        ExtractIlluminaBarcodes \
        NUM_PROCESSORS=16 \
        OUTPUT_DIR=${PICARD_BARCODE_DIR} \
        TMP_DIR=${TMP_FOLDER} \
        BASECALLS_DIR=${FLOWCELL_DIRECTORY}/Data/Intensities/BaseCalls \
        LANE=001 \
        READ_STRUCTURE=36T9B9B36T \
        BARCODE_FILE=${BARCODE_FILE_PATH} \
        METRICS_FILE=${OUTPUT_FOLDER}/metric.txt \
        2> ${OUTPUT_FOLDER}/ExtractIlluminaBarcodes.log; 
    } 2>> $LOG_FILE
};
function run_picard_IlluminaBasecallsToFastq_flat() {
    OUTPUT_FOLDER="${BASE_FOLDER}/basecall/HK5NHBGXX_1"

    rm -rf "$OUTPUT_FOLDER"
    mkdir -p "$OUTPUT_FOLDER"

    echo -e "\nbenchmark IlluminaBasecallsToFastq basecall" >> $LOG_FILE
    clear_os_cache

    {
        time \
        /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx4g \
        -jar /opt/local/tubo/picard-tools-2.7.1/picard.jar \
        IlluminaBasecallsToFastq \
        TMP_DIR=${TMP_FOLDER} \
        NUM_PROCESSORS=16 \
        OUTPUT_PREFIX=${OUTPUT_FOLDER} \
        APPLY_EAMSS_FILTER=false \
        READ_NAME_FORMAT=CASAVA_1_8 \
        FLOWCELL_BARCODE=${FLOWCELL_ID} \
        MAX_READS_IN_RAM_PER_TILE=100000 \
        COMPRESS_OUTPUTS=true \
        BASECALLS_DIR=${FLOWCELL_DIRECTORY}/Data/Intensities/BaseCalls \
        LANE=001 \
        INCLUDE_NON_PF_READS=false \
        READ_STRUCTURE=36T9B9B36T \
        MACHINE_NAME=NB501157 \
        FORCE_GC=true \
        RUN_BARCODE=4 \
        MINIMUM_QUALITY=2 \
        > ${OUTPUT_FOLDER}/IlluminaBasecallsToFastq.log; 
    } 2>> $LOG_FILE
};
function run_picard_IlluminaBasecallsToFastq_demultiplex() {
    OUTPUT_FOLDER="${BASE_FOLDER}/fastq_demultiplexed"
    MULTIPLEX_PARAMS_PATH="${CONFIG_DIRECTORY}/picard_multiplex_params.tsv"

    rm -rf "$OUTPUT_FOLDER"
    mkdir -p "$OUTPUT_FOLDER"

    clear_os_cache

    echo -e "\nbenchmark IlluminaBasecallsToFastq demultiplexed" >> $LOG_FILE
    { 
        time \
        /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx4g \
        -jar /opt/local/tubo/picard-tools-2.7.1/picard.jar \
        IlluminaBasecallsToFastq \
        TMP_DIR=${TMP_FOLDER} \
        BARCODES_DIR=${PICARD_BARCODE_DIR} \
        APPLY_EAMSS_FILTER=false \
        READ_NAME_FORMAT=CASAVA_1_8 \
        FLOWCELL_BARCODE=${FLOWCELL_ID} \
        MAX_READS_IN_RAM_PER_TILE=100000 \
        COMPRESS_OUTPUTS=true \
        BASECALLS_DIR=${FLOWCELL_DIRECTORY}/Data/Intensities/BaseCalls \
        LANE=001 \
        MULTIPLEX_PARAMS=$MULTIPLEX_PARAMS_PATH \
        INCLUDE_NON_PF_READS=false \
        READ_STRUCTURE=${PICARD_READ_STRUCTURE} \
        MACHINE_NAME=NB501157 \
        FORCE_GC=true \
        RUN_BARCODE=4 \
        IGNORE_UNEXPECTED_BARCODES=true \
        2> ${OUTPUT_FOLDER}/IlluminaBasecallsToFastq.log; 
    } 2>> $LOG_FILE
};
function run_picard_IlluminaBasecallsToSam_demultiplex() {
    OUTPUT_FOLDER="${BASE_FOLDER}/cram_demuxltiplex"
    LIBRARY_PARAMS_PATH="${CONFIG_DIRECTORY}/picard_library_params.tsv"

    rm -rf "$OUTPUT_FOLDER"
    mkdir -p "$OUTPUT_FOLDER"

    clear_os_cache

    echo -e "\nbenchmark IlluminaBasecallsToSam demultiplex" >> $LOG_FILE
    { 
        time \
        /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx4g \
        -jar /opt/local/tubo/picard-tools-2.7.1/picard.jar \
        IlluminaBasecallsToSam \
        TMP_DIR=${TMP_FOLDER} \
        BARCODES_DIR=${PICARD_BARCODE_DIR} \
        APPLY_EAMSS_FILTER=false \
        MAX_READS_IN_RAM_PER_TILE=100000 \
        BASECALLS_DIR=${FLOWCELL_DIRECTORY}/Data/Intensities/BaseCalls \
        LANE=001 \
        LIBRARY_PARAMS=$LIBRARY_PARAMS_PATH \
        INCLUDE_NON_PF_READS=false \
        READ_STRUCTURE=${PICARD_READ_STRUCTURE} \
        FORCE_GC=true \
        RUN_BARCODE=HK5NHBGXX \
        PLATFORM=illumina \
        SEQUENCING_CENTER=CGSB \
        2> ${OUTPUT_FOLDER}/IlluminaBasecallsToSam.log; 
    } 2>> $LOG_FILE
};

rm -rf "$BASE_FOLDER"
mkdir -p "$BASE_FOLDER"
mkdir -p "$PICARD_BARCODE_DIR"
mkdir -p "$TMP_FOLDER"

run_picard_ExtractIlluminaBarcodes
run_picard_IlluminaBasecallsToFastq_flat
# run_picard_IlluminaBasecallsToFastq_demultiplex
# run_picard_IlluminaBasecallsToSam_demultiplex
