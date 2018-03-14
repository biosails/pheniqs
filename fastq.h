/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2017  Lior Galanti
    NYU Center for Genetics and System Biology

    Author: Lior Galanti <lior.galanti@nyu.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PHENIQS_FASTQ_H
#define PHENIQS_FASTQ_H

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <htslib/sam.h>
#include <htslib/cram.h>
#include <htslib/bgzf.h>
#include <htslib/kseq.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/hfile.h>
#include <htslib/thread_pool.h>

#include "error.h"
#include "json.h"
#include "constant.h"
#include "auxiliary.h"
#include "sequence.h"
#include "feed.h"
#include "model.h"

using std::map;
using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::string;
using std::vector;
using std::ostream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::make_pair;
using std::setprecision;

using std::mutex;
using std::recursive_mutex;
using std::condition_variable;
using std::unique_lock;
using std::lock_guard;
using std::thread;

KSEQ_INIT(BGZF*, bgzf_read)

/* FASTQ Record
*/
class FastqRecord {
    FastqRecord(FastqRecord const &) = delete;
    void operator=(FastqRecord const &) = delete;

    public:
        kstring_t sequence;
        kstring_t quality;
        kstring_t name;
        kstring_t comment;
        FastqRecord() :
            sequence({ 0, 0, NULL }),
            quality({ 0, 0, NULL }),
            name({ 0, 0, NULL }),
            comment({ 0, 0, NULL }) {
            ks_terminate(sequence);
            ks_terminate(quality);
            ks_terminate(name);
            ks_terminate(comment);
            clear();
        };
        ~FastqRecord() {
            ks_free(sequence);
            ks_free(quality);
            ks_free(name);
            ks_free(comment);
        };
        inline void decode(const kseq_t* kseq, const uint8_t phred_offset) {
            // read a kseq record and populate the FastqRecord
            clear();

            // copy from kseq_t to record
            kputsn(kseq->name.s, kseq->name.l, &name);
            kputsn(kseq->comment.s, kseq->comment.l, &comment);

            // decode sequence
            ks_resize(&sequence, kseq->seq.l + 1);
            for(size_t i = 0; i < kseq->seq.l; i++) {
                sequence.s[i] = AsciiToAmbiguousBam[uint8_t(kseq->seq.s[i])];
                // kputc(AsciiToAmbiguousBam[uint8_t(kseq->seq.s[i])], &sequence);
            }
            sequence.l = kseq->seq.l;
            sequence.s[sequence.l] = '\0';
            // ks_terminate(sequence);

            // decode quality
            ks_resize(&quality, kseq->qual.l + 1);
            for(size_t i = 0; i < kseq->qual.l; i++) {
                quality.s[i] = kseq->qual.s[i] - phred_offset;
                // kputc(kseq->qual.s[i] - phred_offset, &quality);
            }
            sequence.l = kseq->qual.l;
            sequence.s[sequence.l] = '\0';
            // ks_terminate(sequence);
        };
        inline void decode(const Segment& segment) {
            clear();

            // copy from segment to record
            kputsn((char*)segment.sequence.code, segment.sequence.length, &sequence);
            kputsn((char*)segment.sequence.quality, segment.sequence.length, &quality);
            kputsn(segment.name.s, segment.name.l, &name);
            decode_comment(segment);
        };
        inline void encode(Segment& segment) const {

            // write the FastqRecord to Segment
            segment.sequence.fill((const uint8_t*)(sequence.s), (const uint8_t*)(quality.s), sequence.l);
            kputsn(name.s, name.l, &segment.name);
            kputsn(comment.s, comment.l, &segment.auxiliary.CO);
            segment.auxiliary.FI = 0;
            segment.set_qcfail(false);
            segment.auxiliary.XI = 0;
            switch (segment.platform) {
                case Platform::ILLUMINA: {
                    /*  @HWI-ST911:232:HABDFADXX:1:1101:1224:1932 1:N:0:CGATGT
                        name    0:1:2:3:4:5:6
                        0   Instrument ID   HWI-ST911
                        1   Run count       232
                        2   Flowcell ID     HABDFADXX
                        3   Lane number     1
                        4   Tile number     1101
                        5   x coordinate    1224
                        6   y coordinate    1932

                        comment 7:8:9:10
                        7   Segment number  1               uint8_t
                        8   Filtered        N               N|Y
                        9   control number  0               uint16_t
                        10  Barcode         CGATGT          char*
                    */
                    size_t offset = 0;
                    parse_illumina_segment_index(segment, offset);
                    parse_illumina_filtered(segment, offset);
                    parse_illumina_control(segment, offset);
                    parse_illumina_barcode(segment, offset);
                    break;
                }
                case Platform::CAPILLARY:
                case Platform::LS454:
                case Platform::HELICOS:
                case Platform::ONT:
                case Platform::PACBIO:
                case Platform::SOLID:
                case Platform::IONTORRENT:
                    break;
                default:
                    break;
            };
        };
        inline void encode(kstring_t& buffer, uint8_t& phred_offset) const {
            // encode identifier
            kputc('@', &buffer);
            kputsn(name.s, name.l, &buffer);
            kputc(' ', &buffer);
            kputsn(comment.s, comment.l, &buffer);
            kputc(LINE_BREAK, &buffer);

            // encode sequence
            ks_resize(&buffer, buffer.l + sequence.l + 1);
            for (size_t i = 0; i < sequence.l; i++) {
                buffer.s[buffer.l + i] = BamToAmbiguousAscii[uint8_t(sequence.s[i])];
                // kputc(BamToAmbiguousAscii[uint8_t(sequence.s[i])], &buffer);
            }
            buffer.l += sequence.l;
            kputc(LINE_BREAK, &buffer);

            // encode separator
            kputc('+', &buffer);
            kputc(LINE_BREAK, &buffer);

            // encode quality
            ks_resize(&buffer, buffer.l + quality.l + 1);
            for (size_t i = 0; i < quality.l; i++) {
                buffer.s[buffer.l + i] = quality.s[i] + phred_offset;
                // kputc(quality.s[i] + phred_offset, &buffer);
            }
            buffer.l += quality.l;
            kputc(LINE_BREAK, &buffer);
        };

    private:
        inline void clear() {
            ks_clear(sequence);
            ks_clear(quality);
            ks_clear(name);
            ks_clear(comment);
        };
        inline void decode_comment(const Segment& segment) {
            switch (segment.platform) {
                case Platform::CAPILLARY:
                    break;
                case Platform::LS454:
                    break;
                case Platform::ILLUMINA: {
                    kputuw(segment.auxiliary.FI, &comment);
                    kputc(':', &comment);
                    kputc(segment.get_qcfail()? 'Y' : 'N', &comment);
                    kputc(':', &comment);
                    kputuw(segment.auxiliary.XI, &comment);
                    kputc(':', &comment);
                    kputsn(segment.auxiliary.BC.s, segment.auxiliary.BC.l, &comment);
                    break;
                };
                case Platform::SOLID:
                    break;
                case Platform::HELICOS:
                    break;
                case Platform::IONTORRENT:
                    break;
                case Platform::ONT:
                    break;
                case Platform::PACBIO:
                    break;
                default:
                    break;
            };
        };
        inline void parse_illumina_segment_index(Segment& segment, size_t& offset) const {
            int8_t code = -1;
            size_t value = 0;
            size_t position = offset;
            const char* comment = segment.auxiliary.CO.s;
            while (position < segment.auxiliary.CO.l) {
                char c = *(comment + position);
                position++;
                if (c == ':') {
                    break;
                } else {
                    if (code >= -1) {
                        if (c >= '0' && c <= '9') {
                            value *= 10;
                            value += c - '0';
                            code = 0;
                        } else {
                            code = -2;
                        }
                    }
                }
            }
            if (code < 0) {
                value = 1;
            }
            segment.auxiliary.FI = value;
            offset = position;
        };
        inline void parse_illumina_filtered(Segment& segment, size_t& offset) const {
            int8_t code = -1;
            size_t position = offset;
            const char* comment = segment.auxiliary.CO.s;
            while (position < segment.auxiliary.CO.l) {
                char c = *(comment + position);
                position++;
                if (c == ':') {
                    break;
                } else {
                    if (code == -1) {
                        switch(c) {
                            case 'Y':
                                segment.set_qcfail(true);
                                code = 0;
                                break;
                            case 'N':
                                segment.set_qcfail(false);
                                code = 0;
                                break;
                            default:
                                segment.set_qcfail(false);
                                code = -3;
                                break;
                        }
                    } else {
                        // code equals -2 means there was more than one character to the next separator
                        code = -2;
                    }
                }
            }
            offset = position;
        };
        inline void parse_illumina_control(Segment& segment, size_t& offset) const {
            int8_t code = -1;
            uint16_t value = 0;
            size_t position = offset;
            const char* comment = segment.auxiliary.CO.s;
            while (position < segment.auxiliary.CO.l) {
                char c = *(comment + position);
                position++;
                if (c == ':') {
                    break;

                } else {
                    if (code >= -1) {
                        if (c >= '0' && c <= '9') {
                            value *= 10;
                            value += c - '0';
                            code = 0;
                        } else {
                            code = -2;
                        }
                    }
                }
            }
            if (code < 0) {
                value = 0;
            }
            segment.auxiliary.XI = value;
            offset = position;
        };
        inline void parse_illumina_barcode(Segment& segment, size_t& offset) const {
            size_t position = offset;
            const char* comment = segment.auxiliary.CO.s;
            while (position < segment.auxiliary.CO.l) {
                char c = *(comment + position);
                position++;
                if (c == ' ') {
                    break;
                }
            }
            size_t length = position - offset;
            if(length > 0) {
                kputsn(comment + offset, length, &(segment.auxiliary.BC));
            }
            offset = position;
        };
};
class FastqFeed : public BufferedFeed<FastqRecord> {
    friend class Channel;

    public:
        FastqFeed(FeedSpecification const * specification) :
            BufferedFeed<FastqRecord>(specification),
            bgzf_file(NULL) {
        };
        void open() {
            if(!opened()) {
                switch(direction) {
                    case IoDirection::IN: {
                        bgzf_file = bgzf_hopen(hfile, "r");
                        if(bgzf_file != NULL) {
                            kseq = kseq_init(bgzf_file);
                            bgzf_thread_pool(bgzf_file, thread_pool->pool, thread_pool->qsize);
                        } else {
                            throw IOError("failed to open " + string(url) + " for reading");
                        }
                        break;
                    };
                    case IoDirection::OUT: {
                        if(specification->url.compression() == "gz") {
                            bgzf_file = bgzf_hopen(hfile, "wg");
                        } else {
                            bgzf_file = bgzf_hopen(hfile, "wu");
                        }
                        if(bgzf_file != NULL) {
                            kseq = kseq_init(bgzf_file);
                            bgzf_thread_pool(bgzf_file, thread_pool->pool, thread_pool->qsize);
                        } else {
                            throw IOError("failed to open " + string(url) + " for writing");
                        }
                        break;
                    };
                }
            }
        };
        void close() {
            if(opened()) {
                bgzf_close(bgzf_file);
                bgzf_file = NULL;

                kseq_destroy(kseq);
                kseq = NULL;
            }
        };
        inline bool opened() {
            return bgzf_file != NULL;
        };

    protected:
        BGZF* bgzf_file;
        kseq_t* kseq;
        inline void encode(FastqRecord* record, const Segment& segment) const {
            record->decode(segment);
        };
        inline void decode(const FastqRecord* record, Segment& segment) {
            record->encode(segment);
        };
        inline void fill_buffer() {
            while(buffer->is_not_full()) {
             /* >=0  length of the sequence (normal)
                -1   end-of-file
                -2   truncated quality string */
                if(kseq_read(kseq) < 0) {
                    end_of_file = true;
                    break;
                } else {
                    buffer->vacant()->decode(kseq, phred_offset);
                    buffer->increment();
                }
            }
        };
        inline void empty_buffer() {
            /*  encode all fastq records in the buffer to
                a string buffer and write them together to the stream */
            ks_clear(kbuffer);
            while(buffer->is_not_empty()) {
                FastqRecord* record = buffer->next();
                record->encode(kbuffer, phred_offset);
                buffer->decrement();
            }

            if(bgzf_write(bgzf_file, kbuffer.s, kbuffer.l) < 0) {
                throw IOError("error writing to " + string(url));
            }
        };
};
#endif /* PHENIQS_FASTQ_H */
