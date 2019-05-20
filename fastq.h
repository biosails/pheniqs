/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
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

#include "include.h"
#include "feed.h"

KSEQ_INIT(BGZF*, bgzf_read)

class FastqRecord {
    public:
        FastqRecord(FastqRecord const &) = delete;
        void operator=(FastqRecord const &) = delete;
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
            ks_put_string(kseq->name, name);
            ks_put_string(kseq->comment, comment);

            // decode sequence
            ks_increase_to_size(sequence, kseq->seq.l + 2);
            for(size_t i(0); i < kseq->seq.l; ++i) {
                sequence.s[i] = AsciiToAmbiguousBam[static_cast< uint8_t >(kseq->seq.s[i])];
            }
            sequence.l = kseq->seq.l;
            sequence.s[sequence.l] = '\0';

            // decode quality
            ks_increase_to_size(quality, kseq->qual.l + 2);
            for(size_t i(0); i < kseq->qual.l; ++i) {
                quality.s[i] = kseq->qual.s[i] - phred_offset;
            }
            quality.l = kseq->qual.l;
            quality.s[quality.l] = '\0';
        };
        inline void decode(const Segment& segment) {
            clear();

            // copy from segment to record
            ks_put_string(reinterpret_cast< char* >(segment.code), segment.length, sequence);
            ks_put_string(reinterpret_cast< char* >(segment.quality), segment.length, quality);
            ks_put_string(segment.name.s, segment.name.l, name);
            decode_comment(segment);
        };
        inline void encode(Segment& segment) const {

            // write the FastqRecord to Segment
            segment.fill(reinterpret_cast< uint8_t* >(sequence.s), reinterpret_cast< uint8_t* >(quality.s), static_cast< int32_t >(sequence.l));
            ks_put_string(name, segment.name);
            ks_put_string(comment, segment.auxiliary.CO);
            segment.auxiliary.FI = 0;
            segment.set_qcfail(false);

            #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
            segment.auxiliary.illumina_control_number = 0;
            #endif

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
                    size_t offset(0);
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
        inline void encode(kstring_t& buffer, const uint8_t phred_offset) const {
            // encode identifier
            ks_put_character('@', buffer);
            ks_put_string(name, buffer);
            ks_put_character(' ', buffer);
            ks_put_string(comment, buffer);
            ks_put_character(LINE_BREAK, buffer);

            // encode sequence
            ks_increase_by_size(buffer, sequence.l + 2);
            for(size_t i(0); i < sequence.l; ++i) {
                buffer.s[buffer.l + i] = BamToAmbiguousAscii[static_cast< uint8_t >(sequence.s[i])];
            }
            buffer.l += sequence.l;
            ks_put_character(LINE_BREAK, buffer);

            // encode separator
            ks_put_character('+', buffer);
            ks_put_character(LINE_BREAK, buffer);

            // encode quality
            ks_increase_by_size(buffer, quality.l + 2);
            for(size_t i(0); i < quality.l; ++i) {
                buffer.s[buffer.l + i] = quality.s[i] + phred_offset;
            }
            buffer.l += quality.l;
            ks_put_character(LINE_BREAK, buffer);
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
                    // Sanger sequencing
                    break;
                case Platform::LS454:
                    break;
                case Platform::ILLUMINA: {
                    ks_put_uint32(segment.auxiliary.FI, comment);
                    ks_put_character(':', comment);
                    ks_put_character(segment.qcfail() ? 'Y' : 'N', comment);
                    ks_put_character(':', comment);

                    #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
                    ks_put_uint16(segment.auxiliary.illumina_control_number, comment);
                    #else
                    ks_put_uint16(0, comment);
                    #endif

                    ks_put_character(':', comment);
                    ks_put_string(segment.auxiliary.BC, comment);
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
            int8_t code(-1);
            uint32_t value(0);
            size_t position(offset);
            const char* comment(segment.auxiliary.CO.s);
            while (position < segment.auxiliary.CO.l) {
                char c(*(comment + position));
                ++position;
                if(c == ':') {
                    break;
                } else {
                    if(code >= -1) {
                        if(c >= '0' && c <= '9') {
                            value *= 10;
                            value += c - '0';
                            code = 0;
                        } else {
                            code = -2;
                        }
                    }
                }
            }
            if(code < 0) {
                value = 1;
            }
            segment.auxiliary.FI = value;
            offset = position;
        };
        inline void parse_illumina_filtered(Segment& segment, size_t& offset) const {
            int8_t code(-1);
            size_t position(offset);
            const char* comment(segment.auxiliary.CO.s);
            while (position < segment.auxiliary.CO.l) {
                char c(*(comment + position));
                ++position;
                if(c == ':') {
                    break;
                } else {
                    if(code == -1) {
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
            int8_t code(-1);
            uint16_t value(0);
            size_t position(offset);
            const char* comment(segment.auxiliary.CO.s);
            while (position < segment.auxiliary.CO.l) {
                char c(*(comment + position));
                ++position;
                if(c == ':') {
                    break;

                } else {
                    if(code >= -1) {
                        if(c >= '0' && c <= '9') {
                            value *= 10;
                            value += c - '0';
                            code = 0;
                        } else {
                            code = -2;
                        }
                    }
                }
            }
            if(code < 0) {
                value = 0;
            }
            #if defined(PHENIQS_ILLUMINA_CONTROL_NUMBER)
            segment.auxiliary.illumina_control_number = value;
            #endif
            offset = position;
        };
        inline void parse_illumina_barcode(Segment& segment, size_t& offset) const {
            size_t position(offset);
            const char* comment(segment.auxiliary.CO.s);
            while (position < segment.auxiliary.CO.l) {
                char c(*(comment + position));
                ++position;
                if(c == ' ') {
                    break;
                }
            }
            size_t length(position - offset);
            if(length > 0) {
                ks_put_string(comment + offset, length, segment.auxiliary.BC);
            }
            offset = position;
        };
};

class FastqFeed : public BufferedFeed< FastqRecord > {
    friend class Channel;

    public:
        FastqFeed(const FeedProxy& proxy) :
            BufferedFeed< FastqRecord >(proxy),
            bgzf_file(NULL) {
        };
        void open() override {
            if(!opened()) {
                /*  from htslib bgzf.h
                    mode matching /[rwag][u0-9]+/: 'r' for reading, 'w' for
                    writing, 'a' for appending, 'g' for gzip rather than BGZF
                    compression (with 'w' only), and digit specifies the zlib
                    compression level.
                    Note that there is a distinction between 'u' and '0': the
                    first yields plain uncompressed output whereas the latter
                    outputs uncompressed data wrapped in the zlib format.
                 */
                string mode;
                switch(direction) {
                    case IoDirection::IN: {
                        mode.push_back('r');
                        bgzf_file = bgzf_hopen(hfile, mode.c_str());
                        if(bgzf_file != NULL) {
                            kseq = kseq_init(bgzf_file);
                            // bgzf_thread_pool(bgzf_file, thread_pool->pool, thread_pool->qsize);
                        } else {
                            throw IOError("failed to open " + string(url.path()) + " for reading");
                        }
                        break;
                    };
                    case IoDirection::OUT: {
                        mode.push_back('w');
                        switch(url.compression()) {
                            case FormatCompression::GZIP: {
                                mode.push_back('g');
                                break;
                            };
                            case FormatCompression::BGZF: {
                                break;
                            };
                            case FormatCompression::NONE: {
                                mode.push_back('u');
                                break;
                            };
                            default: {
                                mode.push_back('u');
                                break;
                            };
                        };
                        if(url.compression_level() != CompressionLevel::UNKNOWN) {
                            mode.append(to_string(url.compression_level()));
                        }
                        bgzf_file = bgzf_hopen(hfile, mode.c_str());

                        if(bgzf_file != NULL) {
                            kseq = kseq_init(bgzf_file);
                            bgzf_thread_pool(bgzf_file, thread_pool->pool, thread_pool->qsize);
                        } else {
                            throw IOError("failed to open " + string(url.path()) + " for writing");
                        }
                        break;
                    };
                    default: {
                        break;
                    }
                }
            }
        };
        void close() override {
            if(opened()) {
                bgzf_close(bgzf_file);
                bgzf_file = NULL;

                kseq_destroy(kseq);
                kseq = NULL;
            }
        };
        inline bool opened() override {
            return bgzf_file != NULL;
        };

    protected:
        BGZF* bgzf_file;
        kseq_t* kseq;
        inline void encode(FastqRecord* record, const Segment& segment) const override {
            record->decode(segment);
        };
        inline void decode(const FastqRecord* record, Segment& segment) override {
            record->encode(segment);
        };
        inline void replenish_buffer() override {
            while(opened() && buffer->is_not_full()) {
             /* >=0  length of the sequence (normal)
                -1   end-of-file
                -2   truncated quality string */
                if(kseq_read(kseq) < 0) {
                    close();
                    break;
                } else {
                    buffer->vacant()->decode(kseq, phred_offset);
                    buffer->increment();
                }
            }
        };
        inline void flush_buffer() override {
            /*  encode all fastq records in the buffer to
                a string buffer and write them together to the stream */
            if(buffer->is_not_empty()) {
                ks_clear(kbuffer);
                while(buffer->is_not_empty()) {
                    FastqRecord* record = buffer->next();
                    record->encode(kbuffer, phred_offset);
                    buffer->decrement();
                }
                if(bgzf_write(bgzf_file, kbuffer.s, kbuffer.l) < 0) {
                    throw IOError("error writing to " + string(url.path()));
                }
            }
        };
};
#endif /* PHENIQS_FASTQ_H */
