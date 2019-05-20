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
        int32_t     tid;            chromosome ID, defined by bam_hdr_t
        int32_t     pos;            0-based leftmost coordinate
        uint16_t    bin;            bin calculated by bam_reg2bin()
        uint8_t     qual;           mapping quality
        uint8_t     l_qname;        length of the query name
        uint16_t    flag;           bitwise flag
        uint8_t     unused1;
        uint8_t     l_extranul;     length of extra NULs between qname & cigar (for alignment)
        uint32_t    n_cigar;        number of CIGAR operations
        int32_t     l_qseq;         length of the query sequence (read)
        int32_t     mtid;           chromosome ID of next read in template, defined by bam_hdr_t
        int32_t     mpos;           0-based leftmost coordinate of next read in template
        int32_t     isize;
    } bam1_core_t;

    typedef struct {
        bam1_core_t core;
        int         l_data;
        uint32_t    m_data;
        uint8_t*    data;
        uint64_t    id;
    } bam1_t;

    typedef struct {
        int32_t n_targets;          number of reference sequences
        int32_t ignore_sam_err;
        uint32_t l_text;            length of the plain text in the header
        uint32_t *target_len;       lengths of the reference sequences
        int8_t *cigar_tab;
        char **target_name;         names of the reference sequences
        char *text;                 plain text
        void *sdict;                header dictionary
        uint32_t ref_count;
    } bam_hdr_t;
*/

#define bam1_seq_set_i(s, i, c) ((s)[(i)>>0x1] = ((s)[(i)>>0x1]&0xf<<(((i)&0x1)<<0x2))|(c)<<((~(i)&0x1)<<0x2))

ostream& operator<<(ostream& o, const bam1_t& record);

class HtsFeed : public BufferedFeed< bam1_t > {
    friend class Channel;

    public:
        HtsFeed(const FeedProxy& proxy) :
            BufferedFeed< bam1_t >(proxy),
            head(proxy.head),
            hts_file(NULL),
            hdr(NULL) {
        };
        void open() override {
            if(!opened()) {
                /*  from htslib hts.h
                    mode matching / [rwa][bceguxz0-9]* /
                    With 'r' opens for reading; any further format mode letters are ignored
                    as the format is detected by checking the first few bytes or BGZF blocks
                    of the file.  With 'w' or 'a' opens for writing or appending, with format
                    specifier letters:
                        b  binary format (BAM, BCF, etc) rather than text (SAM, VCF, etc)
                        c  CRAM format
                        g  gzip compressed
                        u  uncompressed
                        z  bgzf compressed
                        [0-9]  zlib compression level

                    and with non-format option letters (for any of 'r'/'w'/'a'):
                        e  close the file on exec(2) (opens with O_CLOEXEC, where supported)
                        x  create the file exclusively (opens with O_EXCL, where supported)
                    Note that there is a distinction between 'u' and '0': the first yields
                    plain uncompressed output whereas the latter outputs uncompressed data
                    wrapped in the zlib format.

                    [rw]b  .. compressed BCF, BAM, FAI
                    [rw]bu .. uncompressed BCF
                    [rw]z  .. compressed VCF
                    [rw]   .. uncompressed VCF
                */
                string mode;
                switch(direction) {
                    case IoDirection::IN: {
                        mode.push_back('r');
                        hts_file = hts_hopen(hfile, url.hfile_name(), mode.c_str());
                        if(hts_file != NULL) {
                            hts_set_thread_pool(hts_file, thread_pool);
                            hdr = sam_hdr_read(hts_file);
                            if(hdr != NULL) {
                                head.decode(hdr);
                            } else {
                                throw IOError("failed to read hts header");
                            }
                        } else { throw IOError("failed to open hfile " + string(url.path()) + " for reading"); }
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
                                mode.push_back('z');
                                break;
                            };
                            case FormatCompression::NONE: {
                                mode.push_back('u');
                                break;
                            };
                            default: {
                                break;
                            };
                        };
                        if(url.compression_level() != CompressionLevel::UNKNOWN) {
                            mode.append(to_string(url.compression_level()));
                        }
                        switch(url.type()) {
                            case FormatType::SAM:
                                hts_file = hts_hopen(hfile, url.hfile_name(), mode.c_str());
                                if(hts_file != NULL) {
                                    hts_file->format.version.major = 1;
                                    hts_file->format.version.minor = 0;
                                }
                                break;
                            case FormatType::BAM:
                                mode.push_back('b');
                                hts_file = hts_hopen(hfile, url.hfile_name(), mode.c_str());
                                if(hts_file != NULL) {
                                    hts_file->format.version.major = 1;
                                    hts_file->format.version.minor = 0;
                                }
                                break;
                            case FormatType::CRAM:
                                mode.push_back('c');
                                hts_file = hts_hopen(hfile, url.hfile_name(), mode.c_str());
                                if(hts_file != NULL) {
                                    hts_file->format.version.major = 3;
                                    hts_file->format.version.minor = 0;
                                }
                                break;
                            default:
                                break;
                        }
                        if(hts_file != NULL) {
                            hdr = bam_hdr_init();
                            hts_set_thread_pool(hts_file, thread_pool);
                            head.set_format_version(hts_file->format);
                            head.encode(hdr);
                            if(sam_hdr_write(hts_file, hdr) < 0) {
                                throw IOError("failed to write SAM header");
                            }
                        } else { throw IOError("failed to open hfile " + string(url.path()) + " for writing"); }
                        break;
                    };
                    default:
                        break;
                }
            }
        };
        void close() override {
            if(opened()) {
                // int hts_close_error(hts_close(hts_file));
                // if(hts_close_error) cerr << hts_close_error << endl;
                hts_close(hts_file);
                hts_file = NULL;

                bam_hdr_destroy(hdr);
                hdr = NULL;
            }
        };
        inline bool opened() override {
            return hts_file != NULL;
        };

    protected:
        HtsHead head;
        htsFile* hts_file;
        bam_hdr_t* hdr;
        inline void encode(bam1_t* record, const Segment& segment) const override {
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
                        bam1_seq_set_i(position, i, segment.code[i]);
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
        inline void decode(const bam1_t* record, Segment& segment) override {
            /*  copy the identifier to the segment
                l_qname is :
                    the number of characters in qname +
                    1 for a terminating \0 +
                    l_extranul which pads it to modulo 4 so is 0 to 3
                    in total we need to copy: l_qname - (l_extranul + 1) characters
            */
            ks_put_string(bam_get_qname(record), record->core.l_qname - record->core.l_extranul - 1, segment.name);

            /* copy the sequence */
            uint8_t* bam_seq(bam_get_seq(record));
            segment.increase_to_size(record->core.l_qseq);
            for(int32_t i(0); i < record->core.l_qseq; ++i) {
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
                for(int32_t i(0); i < record->core.l_qseq; ++i) {
                    reinterpret_cast< uint8_t* >(kbuffer.s)[i] = bam_seqi(position, i);
                }
                kbuffer.l = record->core.l_qseq;
                kbuffer.s[kbuffer.l] = '\0';
                segment.fill(reinterpret_cast< uint8_t* >(kbuffer.s), bam_get_qual(record), record->core.l_qseq);
            */

            segment.flag = record->core.flag;
            segment.auxiliary.decode(record);
        };
        inline void replenish_buffer() override {
            while(opened() && buffer->is_not_full()) {
                if(sam_read1(hts_file, hdr, buffer->vacant()) < 0) {
                    close();
                    break;
                } else {
                    buffer->increment();
                }
            }
        };
        inline void flush_buffer() override {
            while(buffer->is_not_empty()) {
                if(sam_write1(hts_file, hdr, buffer->next()) < 0) {
                    throw IOError("error writing to " + string(url.path()));
                }
                buffer->decrement();
            }
        };
};

#endif /* PHENIQS_HTS_H */
