#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# M00518_0198_000000000-A5KVK_VW_17012014:1:1101:15604:1334	77	*	0	0	*	*	0	0	TTCTTCCTTCCTTCTTCTCTTCCTCCTTTCTTTCCCTCCTTTCCTTTCTTTCTCCTTCTTTCCTCCCCTTTCCTTCCCCCTTTCTTTCTCTCCCTCTCTCCTCCCCCTTCCTCCTCGCCTCCCTCCTCCCCCTCCTCCCCCCCCCTCCCT	9<EDEDCBCHHD>GIC<B.EFHGCI;@EH@I>C0A-DECCDCAII:B<BFIG:GGDE8G:C;@A8G3:IC@IEFFD8<I1I@DHCD5IE<C8IIFI6IGDGE4FH;CDEIIFHIFI8;IE?GIEHHEIIBFB7I(EI/IH4IC3:=94?.	XI:Z:CCCTCCT	YI:Z:@IID@II	XJ:Z:TCGCCGG	YJ:Z:BIFIIEF

import json
import sys
import re

def to_json(node):
    return json.dumps(node, sort_keys=True, ensure_ascii=False, indent=4)

i7_seq_tag = 'XI'
i7_qual_tag = 'XJ'
i5_seq_tag = 'YI'
i5_qual_tag = 'YJ'

sam_tag_ex = re.compile(r'^(?P<TAG>[A-Za-z0-9]{2}):(?P<TYPE>[AifZHB]):(?P<VALUE>.+)$')
sam_record_ex = re.compile (
    '^{}$'.format (
        '\t'.join (
            [
                '(?P<QNAME>[^\t]+)',
                '(?P<FLAG>[^\t]+)',
                '(?P<RNAME>[^\t]+)',
                '(?P<POS>[^\t]+)',
                '(?P<MAPQ>[^\t]+)',
                '(?P<CIGAR>[^\t]+)',
                '(?P<RNEXT>[^\t]+)',
                '(?P<PNEXT>[^\t]+)',
                '(?P<TLEN>[^\t]+)',
                '(?P<SEQ>[^\t]+)',
                '(?P<QUAL>[^\t]+)',
                '(?P<AUX>.*)',
            ]
        )
    )
)

# FASTER : parsing the SAM record as a tuple
# about 40% faster but less clear
def parse_sam_tuple(line):
    record = None
    if line:
        tuple = line.split('\t')
        record = { 'tuple': tuple[0:11] }
        if len(tuple) > 11:
            record['AUX'] = {}
            for tag in tuple[11:]:
                tag_match = sam_tag_ex.search(tag)
                if tag_match:
                    tag_record = tag_match.groupdict()
                    record['AUX'][tag_record['TAG']] = tag_record
    return record

def format_segment_tuple(record):
    return '{}\tTC:i:{}'.format('\t'.join(record['tuple']), str(record['TC']))

def process_sam_tuple():
    try:
        AUX = None
        QNAME = None
        segment_index = None

        for line in sys.stdin:
            record = parse_sam_tuple(line)

            if not QNAME == record['tuple'][0]:
                segment_index = 1
                QNAME = record['tuple'][0]
                AUX = record['AUX']
            else:
                segment_index += 1

            record['TC'] = segment_index
            print(format_segment_tuple(record))

            # if this is the first segmentm, print the two index segments
            if AUX is not None:
                # first turn off 0x40 and 0x80 / first and last segment
                record['tuple'][1] = str(int(record['tuple'][1]) & ~0xc0)

                # print i7
                segment_index += 1
                record['tuple'][9] = AUX[i7_seq_tag]['VALUE']
                record['tuple'][10] = AUX[i7_qual_tag]['VALUE']
                record['TC'] = segment_index
                print(format_segment_tuple(record))

                # print i5
                segment_index += 1
                record['tuple'][9] = AUX[i5_seq_tag]['VALUE']
                record['tuple'][10] = AUX[i5_qual_tag]['VALUE']
                record['TC'] = segment_index
                print(format_segment_tuple(record))

                AUX = None

            # print(to_json(record))

    except json.decoder.JSONDecodeError as e:
        print(e)
        sys.exit(1)

    sys.exit(0)

# parsing the SAM record as a dictionary
def parse_sam_record(line):
    record = None
    sam_match = sam_record_ex.search(line)
    if sam_match:
        record = sam_match.groupdict()
        if record['AUX']:
            AUX = {}
            for tag in record['AUX'].split('\t'):
                tag_match = sam_tag_ex.search(tag)
                if tag_match:
                    tag_record = tag_match.groupdict()
                    AUX[tag_record['TAG']] = tag_record
        del record['AUX']
        if AUX:
            record['AUX'] = AUX
    return record

def format_segment(record):
    return '{QNAME}\t{FLAG}\t{RNAME}\t{POS}\t{MAPQ}\t{CIGAR}\t{RNEXT}\t{PNEXT}\t{TLEN}\t{SEQ}\t{QUAL}\tTC:i:{TC}'.format(**record)

def process_sam_input():
    try:
        AUX = None
        QNAME = None
        segment_index = None

        for line in sys.stdin:
            record = parse_sam_record(line)

            if not QNAME == record['QNAME']:
                segment_index = 1
                QNAME = record['QNAME']
                AUX = record['AUX']
            else:
                segment_index += 1

            record['TC'] = segment_index
            print(format_segment(record))

            # if this is the first segmentm, print the two index segments
            if AUX is not None:
                # first turn off 0x40 and 0x80 / first and last segment
                record['FLAG'] = str(int(record['FLAG']) & ~0xc0)

                # print i7
                segment_index += 1
                record['SEQ'] = AUX[i7_seq_tag]['VALUE']
                record['QUAL'] = AUX[i7_qual_tag]['VALUE']
                record['TC'] = segment_index
                print(format_segment(record))

                # print i5
                segment_index += 1
                record['SEQ'] = AUX[i5_seq_tag]['VALUE']
                record['QUAL'] = AUX[i5_qual_tag]['VALUE']
                record['TC'] = segment_index
                print(format_segment(record))

                AUX = None

            # print(to_json(record))

    except json.decoder.JSONDecodeError as e:
        print(e)
        sys.exit(1)

    sys.exit(0)

process_sam_tuple()
