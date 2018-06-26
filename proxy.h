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

#ifndef PHENIQS_PROXY_H
#define PHENIQS_PROXY_H

#include "include.h"

#include <zlib.h>
#include <htslib/hfile.h>
#include "url.h"
#include "atom.h"

const ssize_t PEEK_BUFFER_CAPACITY(4096);
const int DEFAULT_FEED_CAPACITY(60);
const int DEFAULT_FEED_RESOLUTION(60);

enum class FormatKind : uint8_t {
    UNKNOWN,
    DEV_NULL,
    FASTQ,
    HTS,
};
void to_string(const FormatKind& value, string& result);
bool from_string(const char* value, FormatKind& result);
void to_kstring(const FormatKind& value, kstring_t& result);
bool from_string(const string& value, FormatKind& result);
ostream& operator<<(ostream& o, const FormatKind& value);
void encode_key_value(const string& key, const FormatKind& value, Value& container, Document& document);


class FeedProxy {
    friend ostream& operator<<(ostream& o, const FeedProxy& proxy);

    public:
        int32_t index;
        URL url;
        IoDirection direction;
        uint8_t phred_offset;
        hFILE* hfile;
        int capacity;
        int resolution;
        Platform platform;
        unordered_map< string, const HeadPGAtom > program_by_id;
        unordered_map< string, const HeadRGAtom > read_group_by_id;
        FeedProxy(const Value& ontology) :
            index(decode_value_by_key< int32_t >("index", ontology)),
            url(decode_value_by_key< URL >("url", ontology)),
            direction(decode_value_by_key< IoDirection >("direction", ontology)),
            phred_offset(decode_value_by_key< uint8_t >("phred offset", ontology)),
            hfile(NULL),
            capacity(decode_value_by_key< int32_t >("capacity", ontology)),
            resolution(decode_value_by_key< int32_t >("resolution", ontology)),
            platform(decode_value_by_key< Platform >("platform", ontology)) {
        };
        inline bool is_dev_null() const {
            return url.is_dev_null();
        };
        inline bool is_stdin() const {
            return url.is_stdin();
        };
        inline bool is_stdout() const {
            return url.is_stdout();
        };
        inline bool is_stderr() const {
            return url.is_stderr();
        };
        inline FormatKind kind() const {
            if(is_dev_null()) {
                return FormatKind::DEV_NULL;
            } else {
                switch(url.type()) {
                    case FormatType::SAM:
                    case FormatType::BAM:
                    case FormatType::BAI:
                    case FormatType::CRAM:
                    case FormatType::CRAI:
                    case FormatType::VCF:
                    case FormatType::BCF:
                    case FormatType::CSI:
                    case FormatType::GZI:
                    case FormatType::TBI:
                    case FormatType::BED:
                        return FormatKind::HTS;
                        break;
                    case FormatType::FASTQ:
                        return FormatKind::FASTQ;
                        break;
                    default:
                        return FormatKind::UNKNOWN;
                        break;
                }
            }
        };
        void set_capacity(const int& capacity) {
            if(capacity != this->capacity) {
                int aligned(static_cast< int >(capacity / resolution) * resolution);
                if(aligned < capacity) {
                    aligned += resolution;
                }
                this->capacity = aligned;
            }
        };
        void set_resolution(const int& resolution) {
            if(resolution != this->resolution) {
                int aligned(static_cast< int >(capacity / resolution) * resolution);
                if(aligned < capacity) {
                    aligned += resolution;
                }
                this->resolution = resolution;
                this->capacity = aligned;
            }
        };
        void register_rg(const HeadRGAtom& rg) {
            string key(rg);
            auto record = read_group_by_id.find(key);
            if(record == read_group_by_id.end()) {
                read_group_by_id.emplace(make_pair(key, HeadRGAtom(rg)));
            }
        };
        void register_pg(const HeadPGAtom& pg) {
            string key(pg);
            auto record = program_by_id.find(key);
            if(record == program_by_id.end()) {
                program_by_id.emplace(make_pair(key, HeadPGAtom(pg)));
            }
        };
        void probe() {
            if(!url.is_dev_null()) {
                switch(direction) {
                    case IoDirection::IN: {
                        /*  Probe input file

                            Here you can potentially use hfile to probe the file
                            and verify file format and potentially examine the first read
                        */
                        hfile = hopen(url.c_str(), "r");
                        if(url.type() == FormatType::UNKNOWN) {
                            ssize_t peeked(0);
                            unsigned char* buffer(NULL);
                            const ssize_t buffer_capacity(PEEK_BUFFER_CAPACITY);
                            if((buffer = static_cast< unsigned char* >(malloc(buffer_capacity))) == NULL) {
                                throw OutOfMemoryError();
                            }
                            htsFormat format;
                            if(!hts_detect_format(hfile, &format)) {
                                switch (format.format) {
                                    case htsExactFormat::sam:
                                        url.set_type(FormatType::SAM);
                                        break;
                                    case htsExactFormat::bam:
                                        url.set_type(FormatType::BAM);
                                        break;
                                    case htsExactFormat::bai:
                                        url.set_type(FormatType::BAI);
                                        break;
                                    case htsExactFormat::cram:
                                        url.set_type(FormatType::CRAM);
                                        break;
                                    case htsExactFormat::crai:
                                        url.set_type(FormatType::CRAI);
                                        break;
                                    case htsExactFormat::vcf:
                                        url.set_type(FormatType::VCF);
                                        break;
                                    case htsExactFormat::bcf:
                                        url.set_type(FormatType::BCF);
                                        break;
                                    case htsExactFormat::csi:
                                        url.set_type(FormatType::CSI);
                                        break;
                                    case htsExactFormat::gzi:
                                        url.set_type(FormatType::GZI);
                                        break;
                                    case htsExactFormat::tbi:
                                        url.set_type(FormatType::TBI);
                                        break;
                                    case htsExactFormat::bed:
                                        url.set_type(FormatType::BED);
                                        break;
                                    default:
                                        url.set_type(FormatType::UNKNOWN);
                                        break;
                                }
                            }
                            if(url.type() == FormatType::SAM) {
                                peeked = hpeek(hfile, buffer, buffer_capacity);
                                if(peeked > 0) {
                                    switch (format.compression) {
                                        case htsCompression::gzip:
                                        case htsCompression::bgzf: {
                                            unsigned char* decompressed_buffer(NULL);
                                            if((decompressed_buffer = static_cast< unsigned char* >(malloc(buffer_capacity))) == NULL) {
                                                throw OutOfMemoryError();
                                            }
                                            z_stream zstream;
                                            zstream.zalloc = NULL;
                                            zstream.zfree = NULL;
                                            zstream.next_in = buffer;
                                            zstream.avail_in = static_cast< unsigned >(peeked);
                                            zstream.next_out = decompressed_buffer;
                                            zstream.avail_out = buffer_capacity;
                                            if(inflateInit2(&zstream, 31) == Z_OK) {
                                                while(zstream.total_out < buffer_capacity) {
                                                    if(inflate(&zstream, Z_SYNC_FLUSH) != Z_OK) break;
                                                }
                                                inflateEnd(&zstream);
                                                memcpy(buffer, decompressed_buffer, zstream.total_out);
                                                peeked = zstream.total_out;
                                            } else {
                                                peeked = 0;
                                            }
                                            free(decompressed_buffer);
                                            break;
                                        };
                                        case htsCompression::no_compression:
                                            break;
                                        default:
                                            throw InternalError("unknown compression");
                                            break;
                                    }
                                }
                                if(peeked > 0) {
                                    size_t state(0);
                                    char* position(reinterpret_cast< char * >(buffer));
                                    char* end(position + peeked);
                                    while(position < end && position != NULL) {
                                        if(state == 0) {
                                            if(*position == '\n') {
                                                ++position;
                                            } else {
                                                if(*position == '@') {
                                                    state = 1;
                                                } else {
                                                    break;
                                                }
                                            }
                                        } else if(state == 1) {
                                            if((*position >= 'A' && *position <= 'Z') || (*position >= 'a' && *position <= 'z')) {
                                                state = 2;
                                            } else {
                                                break;
                                            }
                                        } else if(state == 2) {
                                            if(*position == '+' && position < end && *(position + 1) == '\n') {
                                                url.set_type(FormatType::FASTQ);
                                            }
                                            break;
                                        }
                                        if((position = strchr(position, '\n')) != NULL) ++position;
                                    }
                                }
                            }
                            free(buffer);
                        }
                        break;
                    };
                    case IoDirection::OUT: {
                        hfile = hopen(url.c_str(), "w");
                        break;
                    };
                    default:
                        break;
                }
            }
        };
};
ostream& operator<<(ostream& o, const FeedProxy& proxy);
template<> FeedProxy decode_value< FeedProxy >(const Value& container);
template<> list< FeedProxy > decode_value_by_key< list< FeedProxy > >(const Value::Ch* key, const Value& container);

#endif /* PHENIQS_PROXY_H */
