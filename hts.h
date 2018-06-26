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

#ifndef PHENIQS_HTS_H
#define PHENIQS_HTS_H

#include "include.h"

#include "feed.h"

/*
    typedef struct {
        int32_t     tid;
        int32_t     pos;
        uint16_t    bin;
        uint8_t     qual;
        uint8_t     l_qname;
        uint16_t    flag;
        uint8_t     unused1;
        uint8_t     l_extranul;
        uint32_t    n_cigar;    *
        int32_t     l_qseq;
        int32_t     mtid;
        int32_t     mpos;
        int32_t     isize;
    } bam1_core_t;

    typedef struct {
        bam1_core_t core;
        int         l_data;
        uint32_t    m_data;
        uint8_t*    data;
        uint64_t    id;
    } bam1_t;
*/

#define bam1_seq_seti(s, i, c) ((s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2))

/*  HTS header */
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

ostream& operator<<(ostream& o, const bam1_t& value);

class HtsFeed : public BufferedFeed< bam1_t > {
    friend class Channel;

    public:
        HtsFeed(const FeedProxy& proxy) :
            BufferedFeed< bam1_t >(proxy),
            hts_file(NULL) {

            header.hd.set_alignment_sort_order(HtsSortOrder::UNKNOWN);
            header.hd.set_alignment_grouping(HtsGrouping::QUERY);
            for(const auto& record : proxy.program_by_id) {
                header.add_program(record.second);
            }
            for(const auto& record : proxy.read_group_by_id) {
                header.add_read_group(record.second);
            }
        };
        void open() {
            if(!opened()) {
                switch(direction) {
                    case IoDirection::IN: {
                        hts_file = hts_hopen(hfile, url.c_str(), "r");
                        if(hts_file != NULL) {
                            hts_set_thread_pool(hts_file, thread_pool);
                            header.decode(hts_file);
                        } else {
                            throw IOError("failed to open " + string(url) + " for reading");
                        }
                        break;
                    };
                    case IoDirection::OUT: {
                        switch(url.type()) {
                            case FormatType::SAM:
                                hts_file = hts_hopen(hfile, url.c_str(), "w");
                                if(hts_file) {
                                    hts_file->format.version.major = 1;
                                    hts_file->format.version.minor = 0;
                                }
                                break;
                            case FormatType::BAM:
                                hts_file = hts_hopen(hfile, url.c_str(), "wb");
                                if(hts_file) {
                                    hts_file->format.version.major = 1;
                                    hts_file->format.version.minor = 0;
                                }
                                break;
                            case FormatType::CRAM:
                                hts_file = hts_hopen(hfile, url.c_str(), "wc");
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
                            throw IOError("failed to open " + string(url) + " for writing");
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
            /*
                The total size of a bam1_t record is an int32_t
                bam1_t.l_data =
                bam1_t.core.l_qname +               // uint8_t : length of the \0 terminated qname, \0 padded to modulo 4
                ((bam1_t.core.l_qseq + 1) >> 1) +   // nucleotide sequence in 4 bit BAM encoding
                bam1_t.core.l_qseq +                // quality sequence in ASCII
                (bam1_t.core.n_cigar << 2) +        // 32 bit per cigar operation
                bam1_t.l_aux                        // auxiliary tags added later
            */
            int32_t i;
            uint32_t l_data;
            int32_t qname_nuls(4 - segment.name.l % 4);
            record->core.l_qname = segment.name.l + qname_nuls;
            if(record->core.l_qname <= numeric_limits< uint8_t >::max()) {
                record->core.flag = segment.flag;
                record->core.l_extranul = static_cast< uint8_t >(qname_nuls - 1);
                record->core.l_qseq = static_cast< int32_t >(segment.length);
                l_data = record->core.l_qname + (record->core.n_cigar << 2) + ((segment.length + 1) >> 1) + segment.length;
                if(l_data <= numeric_limits< int32_t >::max()) {
                    record->l_data = static_cast< int32_t >(l_data);
                    if(record->m_data < l_data) {
                        record->m_data = l_data;
                        kroundup32(record->m_data);
                        if((record->data = static_cast< uint8_t* >(realloc(record->data, record->m_data))) == NULL) {
                            throw OutOfMemoryError();
                        }
                    }

                    uint8_t* position(record->data);

                    // write identifier
                    memcpy(position, segment.name.s, segment.name.l);
                    position += segment.name.l;
                    for(i = 0; i < qname_nuls; ++i) {
                        *position = '\0';
                        ++position;
                    }

                    // write cigar string
                    if(record->core.n_cigar > 0) {
                        // memcpy(position, cigar, record->core.n_cigar * 4);
                        position += (record->core.n_cigar << 2);
                    }

                    // encode nucleotide byte BAM numeric encoding into nybble BAM numeric encoding
                    for(i = 0; i < segment.length; ++i) {
                        bam1_seq_seti(position, i, segment.code[i]);
                    }
                    position += ((segment.length + 1) >> 1);

                    /*  alternative sequence encoding implementation

                        for(i = 0; i + 1 < segment.sequence.length; i += 2) {
                            ++(*position) = (AsciiToAmbiguousBam[segment.sequence.code[i]] << 4) + AsciiToAmbiguousBam[segment.sequence.code[i + 1]];
                        }
                        if(i < segment.sequence.length) {
                            ++(*position) = (AsciiToAmbiguousBam[segment.sequence.code[i]] << 4);
                        }
                    */

                    // encode the quality sequence
                    memcpy(position, segment.quality, segment.length);

                    // encode the auxiliary tag which will update l_data accordingly
                    segment.auxiliary.encode(record);

                } else { throw OverflowError("BAM record must not exceed " + to_string(numeric_limits< int32_t >::max()) + " bytes"); }
            } else { throw OverflowError("qname must not exceed 255 characters"); }
        };
        inline void decode(const bam1_t* record, Segment& segment) {
            /* copy the identifier to the segment */
            ks_put_string(bam_get_qname(record), record->core.l_qname, segment.name);

            /* copy the sequence */
            uint8_t* bam_seq(bam_get_seq(record));
            segment.increase_to_size(record->core.l_qseq);
            for(int32_t i = 0; i < record->core.l_qseq; ++i) {
                /* pad 4bit BAM numeric encoding to 8bit */
                segment.code[i] = bam_seqi(bam_seq, i);
            }
            segment.code[record->core.l_qseq] = '\0';

            /* copy the quality */
            memcpy(segment.quality, bam_get_qual(record), record->core.l_qseq);
            segment.quality[record->core.l_qseq] = '\0';

            /* assign the sequence length */
            segment.length = record->core.l_qseq;

            /*  copy sequence using a buffer

                ks_clear(kbuffer);
                ks_increase_to_size(kbuffer, record->core.l_qseq + 2);
                uint8_t* position(bam_get_seq(record));
                for(int32_t i = 0; i < record->core.l_qseq; ++i) {
                    reinterpret_cast< uint8_t* >(kbuffer.s)[i] = bam_seqi(position, i);
                }
                kbuffer.l = record->core.l_qseq;
                kbuffer.s[kbuffer.l] = '\0';
                segment.fill(reinterpret_cast< uint8_t* >(kbuffer.s), bam_get_qual(record), record->core.l_qseq);
            */

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
                    throw IOError("error writing to " + string(url));
                }
                buffer->decrement();
            }
        };
};
#endif /* PHENIQS_HTS_H */
