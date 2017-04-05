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


FLOWCELL="$1"
PHIX_SAMPLE_PATH="HK5NHBGXX_phix"
QUALITY_SAMPLE_PATH="HK5NHBGXX_barcode_quality"
PHIX_QUALITY_SAMPLE_PATH="HK5NHBGXX_phix_barcode_quality"
PHENIQS_CONFIG="${FLOWCELL}_demultiplex.json"
REPORT_PATH="${FLOWCELL}_measure.json"
LOG_DIR="log"
BASE_REMOTE_URL="http://albireo.bio.nyu.edu/~lg/secret/f5841601-b547-4a9c-aa02-70075c9b120a"

DEMULTIPLEXED_CRAM="${FLOWCELL}.cram"

function unpack() {
    if [[ ! -f "${PHIX_SAMPLE_PATH}" ]]; then
        if [[ ! -f "${PHIX_SAMPLE_PATH}.bz2" ]]; then
            curl -O "${BASE_REMOTE_URL}/${PHIX_SAMPLE_PATH}.bz2"
        fi

        if [[ -f "${PHIX_SAMPLE_PATH}.bz2" ]]; then
            echo "uncompressing ${PHIX_SAMPLE_PATH}.bz2 file";
            bzip2 -dk "${PHIX_SAMPLE_PATH}.bz2"
        else
            echo "missing ${PHIX_SAMPLE_PATH} file";
            exit(1)
        fi
    fi

    if [[ ! -f "${QUALITY_SAMPLE_PATH}" ]]; then
        if [[ ! -f "${QUALITY_SAMPLE_PATH}.bz2" ]]; then
            curl -O "$BASE_REMOTE_URL/${QUALITY_SAMPLE_PATH}.bz2"
        fi

        if [[ -f "${QUALITY_SAMPLE_PATH}.bz2" ]]; then
            echo "uncompressing ${QUALITY_SAMPLE_PATH}.bz2 file";
            bzip2 -dk "${QUALITY_SAMPLE_PATH}.bz2"
        else
            echo "missing ${QUALITY_SAMPLE_PATH} file";
            exit(1)
        fi
    fi

    if [[ ! -f "${PHIX_QUALITY_SAMPLE_PATH}" ]]; then
        if [[ ! -f "${PHIX_QUALITY_SAMPLE_PATH}.bz2" ]]; then
            curl -O "$BASE_REMOTE_URL/${PHIX_QUALITY_SAMPLE_PATH}.bz2"
        fi

        if [[ -f "${PHIX_QUALITY_SAMPLE_PATH}.bz2" ]]; then
            echo "uncompressing ${PHIX_QUALITY_SAMPLE_PATH}.bz2 file";
            bzip2 -dk "${PHIX_QUALITY_SAMPLE_PATH}.bz2"
        else
            echo "missing ${PHIX_QUALITY_SAMPLE_PATH} file";
            exit(1)
        fi
    fi
}

function clean_fastq() {
    echo "removing simulated FASTQ" && rm -f *.fastq
}

function clean_report() {
    echo "removing report" && rm -f "${REPORT_PATH}"
}

function clean_log() {
    [[ -d "${LOG_DIR}" ]] && echo "removing ${LOG_DIR}" && rm -rf "${LOG_DIR}"
}

function clean() {
    clean_fastq
    clean_report
    clean_log
}

function analyze() {
    NATURE="$1"
    MODEL_PATH="${FLOWCELL}_${NATURE}_model.json"
    SIMULATED_FASTQ="${FLOWCELL}_${NATURE}.fastq"

    # create a simulated fastq
    if [[ ! -f "${SIMULATED_FASTQ}" ]]; then
        echo "simulating FASTQ ${SIMULATED_FASTQ}";
        time ./simulate.py simulate "${MODEL_PATH}" |pv > "${SIMULATED_FASTQ}"
    fi

    for CONFIDENCE in 85 90 95 98 99 995 998 999 9999 99999; do
        PHENIQS_LOG="${LOG_DIR}/demux_${NATURE}_${CONFIDENCE}.json"

        [ -f "${DEMULTIPLEXED_CRAM}" ] && rm -f "${DEMULTIPLEXED_CRAM}";
        [ -f "${PHENIQS_LOG}" ] && rm -f "${PHENIQS_LOG}";
        
        [ -f "${PHENIQS_CONFIG}" ] && \
        time cat "${SIMULATED_FASTQ}" | pheniqs demux --config "${PHENIQS_CONFIG}" --confidence "0.${CONFIDENCE}" > "${PHENIQS_LOG}"
        
        [ -f "${DEMULTIPLEXED_CRAM}" ] && time samtools view "${DEMULTIPLEXED_CRAM}" |\
        ./simulate.py analyze "${MODEL_PATH}" "${REPORT_PATH}" "0.${CONFIDENCE}";
    done
}

if [[ $FLOWCELL ]]; then
    unpack
    clean_fastq
    clean_report
    clean_log

    [[ ! -d "$LOG_DIR" ]] && echo "creating log directory ${LOG_DIR}" && mkdir -p "${LOG_DIR}"

    analyze "natural"
else
    echo "must specify a flowcell";
fi
