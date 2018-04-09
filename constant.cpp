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

#include "constant.h"

/*  Environment constants */

void to_string(const FormatKind& value, string& result) {
    switch(value) {
        case FormatKind::FASTQ: result.assign("FASTQ");      break;
        case FormatKind::HTS:   result.assign("HTS");        break;
        default:                result.assign("UNKNOWN");    break;
    }
};
bool from_string(const char* value, FormatKind& result) {
         if(value == NULL)              result = FormatKind::UNKNOWN;
    else if(!strcmp(value, "FASTQ"))    result = FormatKind::FASTQ;
    else if(!strcmp(value, "HTS"))      result = FormatKind::HTS;
    else                                result = FormatKind::UNKNOWN;

    return (result == FormatKind::UNKNOWN ? false : true);
};
void to_kstring(const FormatKind& value, kstring_t& result) {
    ks_clear(result);
    string string_value;
    to_string(value, string_value);
    kputs(string_value.c_str(), &result);
};
bool from_string(const string& value, FormatKind& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const FormatKind& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const FormatKind& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< FormatKind >(const Value::Ch* key, FormatKind& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};

void to_string(const Decoder& value, string& result) {
    switch(value) {
        case Decoder::UNKNOWN:      result.assign("unknown");    break;
        case Decoder::MDD:          result.assign("mdd");        break;
        case Decoder::PAMLD:        result.assign("pamld");      break;
        case Decoder::BENCHMARK:    result.assign("benchmark");  break;
        default:                                                 break;
    }
};
bool from_string(const char* value, Decoder& result) {
         if(value == NULL)                  result = Decoder::UNKNOWN;
    else if(!strcmp(value, "mdd"))          result = Decoder::MDD;
    else if(!strcmp(value, "pamld"))        result = Decoder::PAMLD;
    else if(!strcmp(value, "benchmark"))    result = Decoder::BENCHMARK;
    else                                    result = Decoder::UNKNOWN;

    return (result == Decoder::UNKNOWN ? false : true);
};
void to_kstring(const Decoder& value, kstring_t& result) {
    ks_clear(result);
    string string_value;
    to_string(value, string_value);
    kputs(string_value.c_str(), &result);
};
bool from_string(const string& value, Decoder& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const Decoder& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
bool encode_key_value(const string& key, const Decoder& value, Value& container, Document& document) {
    if(value != Decoder::UNKNOWN) {
        string string_value;
        to_string(value, string_value);
        Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.RemoveMember(key.c_str());
        container.AddMember(k.Move(), v.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
template<> bool decode_value_by_key< Decoder >(const Value::Ch* key, Decoder& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};

template <> Decoder decode_value_by_key(const Value::Ch* key, const Value& container) {
    Decoder value(Decoder::UNKNOWN);
    decode_value_by_key(key, value, container);
    return value;
};

ostream& operator<<(ostream& o, const LeftTokenOperator& value) {
    switch (value) {
        case LeftTokenOperator::NONE:               o << "none";                break;
        case LeftTokenOperator::REVERSE_COMPLEMENT: o << "reverse complement";  break;
    }
    return o;
};

void to_string(const Platform& value, string& result) {
    switch(value) {
        case Platform::CAPILLARY:   result.assign("CAPILLARY");  break;
        case Platform::LS454:       result.assign("LS454");      break;
        case Platform::ILLUMINA:    result.assign("ILLUMINA");   break;
        case Platform::SOLID:       result.assign("SOLID");      break;
        case Platform::HELICOS:     result.assign("HELICOS");    break;
        case Platform::IONTORRENT:  result.assign("IONTORRENT"); break;
        case Platform::ONT:         result.assign("ONT");        break;
        case Platform::PACBIO:      result.assign("PACBIO");     break;
        default:                    result.assign("UNKNOWN");    break;
    }
};
bool from_string(const char* value, Platform& result) {
         if(value == NULL)                  result = Platform::UNKNOWN;
    else if(!strcmp(value, "CAPILLARY"))    result = Platform::CAPILLARY;
    else if(!strcmp(value, "LS454"))        result = Platform::LS454;
    else if(!strcmp(value, "ILLUMINA"))     result = Platform::ILLUMINA;
    else if(!strcmp(value, "SOLID"))        result = Platform::SOLID;
    else if(!strcmp(value, "HELICOS"))      result = Platform::HELICOS;
    else if(!strcmp(value, "IONTORRENT"))   result = Platform::IONTORRENT;
    else if(!strcmp(value, "ONT"))          result = Platform::ONT;
    else if(!strcmp(value, "PACBIO"))       result = Platform::PACBIO;
    else                                    result = Platform::UNKNOWN;

    return (result == Platform::UNKNOWN ? false : true);
};
void to_kstring(const Platform& value, kstring_t& result) {
    ks_clear(result);
    string string_value;
    to_string(value, string_value);
    kputs(string_value.c_str(), &result);
};
bool from_string(const string& value, Platform& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const Platform& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const Platform& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< Platform >(const Value::Ch* key, Platform& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
template <> Platform decode_value_by_key(const Value::Ch* key, const Value& container) {
    Platform value(Platform::UNKNOWN);
    decode_value_by_key(key, value, container);
    return value;
};

void to_string(const HtsSortOrder& value, string& result) {
    switch(value) {
        case HtsSortOrder::UNSORTED:    result.assign("unsorted");   break;
        case HtsSortOrder::QUERYNAME:   result.assign("queryname");  break;
        case HtsSortOrder::COORDINATE:  result.assign("coordinate"); break;
        default:                        result.assign("unknown");    break;
    }
};
bool from_string(const char* value, HtsSortOrder& result) {
         if(value == NULL)                  result = HtsSortOrder::UNKNOWN;
    else if(!strcmp(value, "unsorted"))     result = HtsSortOrder::UNSORTED;
    else if(!strcmp(value, "queryname"))    result = HtsSortOrder::QUERYNAME;
    else if(!strcmp(value, "coordinate"))   result = HtsSortOrder::COORDINATE;
    else                                    result = HtsSortOrder::UNKNOWN;

    return (result == HtsSortOrder::UNKNOWN ? false : true);
};
void to_kstring(const HtsSortOrder& value, kstring_t& result) {
    ks_clear(result);
    string string_value;
    to_string(value, string_value);
    kputs(string_value.c_str(), &result);
};
bool from_string(const string& value, HtsSortOrder& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const HtsSortOrder& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const HtsSortOrder& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< HtsSortOrder >(const Value::Ch* key, HtsSortOrder& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};

void to_string(const HtsGrouping& value, string& result) {
    switch(value) {
        case HtsGrouping::QUERY:        result.assign("query");      break;
        case HtsGrouping::REFERENCE:    result.assign("reference");  break;
        default:                        result.assign("none");       break;
    }
};
bool from_string(const char* value, HtsGrouping& result) {
         if(value == NULL)                  result = HtsGrouping::NONE;
    else if(!strcmp(value, "query"))        result = HtsGrouping::QUERY;
    else if(!strcmp(value, "reference"))    result = HtsGrouping::REFERENCE;
    else                                    result = HtsGrouping::NONE;

    return (result == HtsGrouping::UNKNOWN ? false : true);
};
void to_kstring(const HtsGrouping& value, kstring_t& result) {
    ks_clear(result);
    string string_value;
    to_string(value, string_value);
    kputs(string_value.c_str(), &result);
};
bool from_string(const string& value, HtsGrouping& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const HtsGrouping& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const HtsGrouping& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< HtsGrouping >(const Value::Ch* key, HtsGrouping& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};

ostream& operator<<(ostream& o, const htsFormatCategory& hts_format_category) {
    switch (hts_format_category) {
        case htsFormatCategory::sequence_data:  o << "sequence data";   break;
        case htsFormatCategory::variant_data:   o << "variant data";    break;
        case htsFormatCategory::index_file:     o << "index file";      break;
        case htsFormatCategory::region_list:    o << "region list";     break;
        default:                                o << "unknown";         break;
    }
    return o;
};

ostream& operator<<(ostream& o, const htsExactFormat& value) {
    switch (value) {
        case htsExactFormat::binary_format: o << "binary format";   break;
        case htsExactFormat::text_format:   o << "text format";     break;
        case htsExactFormat::sam:           o << "sam";             break;
        case htsExactFormat::bam:           o << "bam";             break;
        case htsExactFormat::bai:           o << "bai";             break;
        case htsExactFormat::cram:          o << "cram";            break;
        case htsExactFormat::crai:          o << "crai";            break;
        case htsExactFormat::vcf:           o << "vcf";             break;
        case htsExactFormat::bcf:           o << "bcf";             break;
        case htsExactFormat::csi:           o << "csi";             break;
        case htsExactFormat::gzi:           o << "gzi";             break;
        case htsExactFormat::tbi:           o << "tbi";             break;
        case htsExactFormat::bed:           o << "bed";             break;
        default:                            o << "unknown";         break;
    }
    return o;
};

ostream& operator<<(ostream& o, const htsCompression& value) {
    switch (value) {
        case htsCompression::gzip:      o << "gzip";            break;
        case htsCompression::bgzf:      o << "bgzf";            break;
        case htsCompression::custom:    o << "custom";          break;
        default:                        o << "no compression";  break;
    }
    return o;
};
