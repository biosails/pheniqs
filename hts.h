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

#ifndef PHENIQS_HTS_H
#define PHENIQS_HTS_H

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "error.h"
#include "atom.h"
#include "auxiliary.h"
#include "sequence.h"
#include "feed.h"
#include "specification.h"

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

#define bam1_seq_seti(s, i, c) ((s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2))

/*  HTS header
*/
class HtsHeader {
friend ostream& operator<<(ostream& o, const HtsHeader& head);
HtsHeader(HtsHeader const &) = delete;
void operator=(HtsHeader const &) = delete;

public:
    bam_hdr_t* hdr;
    HeadHDAtom hd;
    vector< HeadCOAtom > comments;
    map< string, const HeadPGAtom > program_by_id;
    map< string, const HeadRGAtom > read_group_by_id;
    HtsHeader();
    ~HtsHeader();
    void assemble();
    void decode(htsFile* hts_file);
    void encode(htsFile* hts_file) const;
    void add_read_group(const HeadRGAtom& read_group);
    void add_program(const HeadPGAtom& program);
    void add_comment(const HeadCOAtom& co);
};
ostream& operator<<(ostream& o, const HtsHeader& header);

class HtsFeed : public BufferedFeed<bam1_t> {
friend class Channel;

public:
    HtsFeed(const FeedSpecification& specification) :
        BufferedFeed<bam1_t>(specification),
        hts_file(NULL) {

        header.hd.set_alignment_sort_order(HtsSortOrder::UNKNOWN);
        header.hd.set_alignment_grouping(HtsGrouping::QUERY);
        for(const auto& record : specification.program_by_id) {
            header.add_program(record.second);
        }
        for(const auto& record : specification.read_group_by_id) {
            header.add_read_group(record.second);
        }
    };
    void open() {
        if(!opened()) {
            switch(specification.direction) {
                case IoDirection::IN: {
                    hts_file = hts_hopen(hfile, specification.url.c_str(), "r");
                    if(hts_file != NULL) {
                        hts_set_thread_pool(hts_file, thread_pool);
                        header.decode(hts_file);
                    } else {
                        throw IOError("failed to open " + string(specification.url) + " for reading");
                    }
                    break;
                };
                case IoDirection::OUT: {
                    switch(specification.url.type()) {
                        case FormatType::SAM:
                            hts_file = hts_hopen(hfile, specification.url.c_str(), "w");
                            if(hts_file) {
                                hts_file->format.version.major = 1;
                                hts_file->format.version.minor = 0;
                            }
                            break;
                        case FormatType::BAM:
                            hts_file = hts_hopen(hfile, specification.url.c_str(), "wb");
                            if(hts_file) {
                                hts_file->format.version.major = 1;
                                hts_file->format.version.minor = 0;
                            }
                            break;
                        case FormatType::CRAM:
                            hts_file = hts_hopen(hfile, specification.url.c_str(), "wc");
                            if(hts_file) {
                                hts_file->format.version.major = 3;
                                hts_file->format.version.minor = 0;
                            }
                            break;
                        default:
                            break;
                    }
                    if(hts_file != NULL) {
                        hts_set_thread_pool(hts_file, thread_pool);
                        header.hd.set_version(&(hts_file->format));
                        header.assemble();
                        header.encode(hts_file);
                    } else {
                        throw IOError("failed to open " + string(specification.url) + " for writing");
                    }
                    break;
                };
                default:
                    break;
            }
        }
    };
    void close() {
        if(opened()) {
            hts_close(hts_file);
            hts_file = NULL;
        }
    };
    inline bool opened() {
        return hts_file != NULL;
    };
    HtsHeader& get_header() {
        return header;
    };

protected:
    HtsHeader header;
    htsFile* hts_file;
    inline void encode(bam1_t* record, const Segment& segment) const {
        record->core.l_qname = segment.name.l + 1;
        record->core.l_qseq = segment.sequence.length;
        record->core.flag = segment.flag;

        uint64_t l_data = record->core.l_qname + ((record->core.l_qseq + 1)>>1) + record->core.l_qseq;
        record->l_data = l_data;

        // increase allocated space if necessary
        if(record->m_data < l_data) {
            record->m_data = l_data;
            kroundup32(record->m_data);
            if((record->data = static_cast< uint8_t* >(realloc(record->data, record->m_data))) == NULL) {
                free(record->data);
                throw InternalError("out of memory");
            }
        }

        // write identifier to data block
        memcpy(bam_get_qname(record), segment.name.s, record->core.l_qname);

        // encode nucleotide byte BAM numeric encoding into nybble BAM numeric encoding
        uint8_t *bam_seq = bam_get_seq(record);
        for(int i = 0; i < record->core.l_qseq; i++) {
            bam1_seq_seti(bam_seq, i, segment.sequence.code[i]);
        }

        // encode the quality sequence
        memcpy(bam_get_qual(record), segment.sequence.quality, record->core.l_qseq);

        segment.auxiliary.encode(record);
    };
    inline void decode(const bam1_t* record, Segment& segment) {
        // decode identifier and write to segment
        ks_put_string(bam_get_qname(record), record->core.l_qname, segment.name);

        // decode nucleotide sequence to buffer
        ks_clear(kbuffer);
        ks_resize(&kbuffer, record->core.l_qseq + 1);

        uint8_t* position = bam_get_seq(record);
        for(int32_t i = 0; i < record->core.l_qseq; i++) {
            // convert from 4bit BAM numeric encoding to 8bit BAM numeric encoding
            kbuffer.s[i] = bam_seqi(position, i);
            // ks_put_character(bam_seqi(position, i), &kbuffer);
        }
        kbuffer.l = record->core.l_qseq;
        kbuffer.s[kbuffer.l] = '\0';
        segment.sequence.fill((uint8_t*)kbuffer.s, bam_get_qual(record), record->core.l_qseq);

        segment.flag = record->core.flag;
        segment.auxiliary.decode(record);
    };
    inline void replenish_buffer() {
        while(opened() && buffer->is_not_full()) {
            if(sam_read1(hts_file, header.hdr, buffer->vacant()) < 0) {
                close();
                break;
            } else {
                buffer->increment();
            }
        }
    };
    inline void flush_buffer() {
        while(buffer->is_not_empty()) {
            if(sam_write1(hts_file, header.hdr, buffer->next()) < 0) {
                throw IOError("error writing to " + string(specification.url));
            }
            buffer->decrement();
        }
    };
};
#endif /* PHENIQS_HTS_H */
