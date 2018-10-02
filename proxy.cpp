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

#include "proxy.h"

FeedProxy::FeedProxy(const Value& ontology) try :
    index(decode_value_by_key< int32_t >("index", ontology)),
    url(decode_value_by_key< URL >("url", ontology)),
    direction(decode_value_by_key< IoDirection >("direction", ontology)),
    phred_offset(decode_value_by_key< uint8_t >("phred offset", ontology)),
    hfile(NULL),
    capacity(decode_value_by_key< int32_t >("capacity", ontology)),
    resolution(decode_value_by_key< int32_t >("resolution", ontology)),
    platform(decode_value_by_key< Platform >("platform", ontology)) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("FeedProxy :: " + error.message);

    } catch(exception& error) {
        throw InternalError("FeedProxy :: " + string(error.what()));
};
void FeedProxy::register_rg(const HeadRGAtom& rg) {
    string key(rg);
    auto record = read_group_by_id.find(key);
    if(record == read_group_by_id.end()) {
        read_group_by_id.emplace(make_pair(key, HeadRGAtom(rg)));
    }
};
void FeedProxy::register_pg(const HeadPGAtom& pg) {
    string key(pg);
    auto record = program_by_id.find(key);
    if(record == program_by_id.end()) {
        program_by_id.emplace(make_pair(key, HeadPGAtom(pg)));
    }
};
void FeedProxy::probe() {
    if(hfile == NULL) {
        if(!is_dev_null()) {
            switch(direction) {
                case IoDirection::IN: {
                    /*  Use hfile to probe the input file and verify format */
                    if(url.is_readable()) {
                        hfile = hopen(url.hfile_name(), "r");
                        if(hfile != NULL) {
                            if(url.type() == FormatType::UNKNOWN) {
                                ssize_t peeked(0);
                                unsigned char* buffer(NULL);
                                const ssize_t buffer_capacity(PEEK_BUFFER_CAPACITY);
                                if((buffer = static_cast< unsigned char* >(malloc(buffer_capacity))) == NULL) {
                                    throw OutOfMemoryError();
                                }

                                /* detect input format with htslib and override the type encoded in the url */
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
                                                z_stream gz_stream;
                                                gz_stream.zalloc = NULL;
                                                gz_stream.zfree = NULL;
                                                gz_stream.next_in = buffer;
                                                gz_stream.avail_in = static_cast< unsigned >(peeked);
                                                gz_stream.next_out = decompressed_buffer;
                                                gz_stream.avail_out = buffer_capacity;
                                                if(inflateInit2(&gz_stream, 31) == Z_OK) {
                                                    while(gz_stream.total_out < buffer_capacity) {
                                                        if(inflate(&gz_stream, Z_SYNC_FLUSH) != Z_OK) break;
                                                    }
                                                    inflateEnd(&gz_stream);
                                                    memcpy(buffer, decompressed_buffer, gz_stream.total_out);
                                                    peeked = gz_stream.total_out;
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
                                                    } else { break; }
                                                }
                                            } else if(state == 1) {
                                                if((*position >= 'A' && *position <= 'Z') || (*position >= 'a' && *position <= 'z')) {
                                                    state = 2;
                                                } else { break; }
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
                        } else { throw IOError("failed to open " + string(url) + " for reading"); }
                    } else { throw IOError("could not open " + string(url) + " for reading"); }
                    break;
                };
                case IoDirection::OUT: {
                    if(url.is_writable()) {
                        hfile = hopen(url.hfile_name(), "w");
                        if(hfile != NULL) {
                            /* more testing */
                        } else { throw IOError("failed to open " + string(url) + " for writing"); }
                    } else { throw IOError("could not open " + string(url) + " for writing"); }
                    break;
                };
                default:
                    break;
            }
        }
    }
};

bool encode_value(const string& key, const FeedProxy& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        Value element(kObjectType);
        encode_key_value("index", value.index, element, document);
        encode_key_value("url", value.url, element, document);
        encode_key_value("direction", value.direction, element, document);
        encode_key_value("platform", value.platform, element, document);
        encode_key_value("capacity", value.capacity, element, document);
        encode_key_value("resolution", value.resolution, element, document);
        encode_key_value("phred offset", value.phred_offset, element, document);
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), element.Move(), document.GetAllocator());
        return true;
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
ostream& operator<<(ostream& o, const FeedProxy& proxy) {
    o << "direction : " << proxy.direction << endl;
    o << "index : " << proxy.index << endl;
    o << "url : " << proxy.url << endl;
    o << "platform : " << proxy.platform << endl;
    o << "capacity : " << proxy.capacity << endl;
    o << "resolution : " << proxy.resolution << endl;
    o << "phred_offset : " << to_string(proxy.phred_offset) << endl;
    o << proxy.url.description();
    return o;
};
template<> FeedProxy decode_value< FeedProxy >(const Value& container) {
    if(container.IsObject()) {
        FeedProxy proxy(container);
        return proxy;
    } else { throw ConfigurationError("feed proxy element must be a dictionary"); }
};
template<> list< FeedProxy > decode_value_by_key< list< FeedProxy > >(const Value::Ch* key, const Value& container) {
    if(container.IsObject()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            if(!reference->value.IsNull()) {
                if(reference->value.IsArray()) {
                    list< FeedProxy > value;
                    if(!reference->value.Empty()) {
                        for(auto& element : reference->value.GetArray()) {
                            value.emplace_back(element);
                        }
                    }
                    return value;
                } else { throw ConfigurationError(string(key) + " is not an array"); }
            } else { throw ConfigurationError(string(key) + " is null"); }
        } else { throw ConfigurationError(string(key) + " not found"); }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
};
