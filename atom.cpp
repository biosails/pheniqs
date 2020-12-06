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

#include "atom.h"

static inline char* copy_until_tag_end(char* source, const char* end, kstring_t& target) {
    ks_clear(target);
    char* position(source);
    if(source < end) {
        while ((*position != '\t' && *position != LINE_BREAK) && position < end) {
            ++position;
        }
        size_t length(position - source);
        if(length > 0) {
            ks_put_string(source, length, target);
        }
    }
    return position;
};
static inline char* copy_until_linebreak(char* source, const char* end, kstring_t& target) {
    ks_clear(target);
    char* position(source);
    if(source < end) {
        while (*position != LINE_BREAK && position < end) {
            ++position;
        }
        size_t length(position - source);
        if(length > 0) {
            ks_put_string(source, length, target);
        }
    }
    return position;
};
static inline char* skip_to_tab(char* source, const char* end) {
    char* position(source);
    while (*position != '\t' && position < end) {
        ++position;
    }
    return position;
};

string to_string(const HtsSortOrder& value) {
    string result;
    switch(value) {
        case HtsSortOrder::UNSORTED:    result.assign("unsorted");   break;
        case HtsSortOrder::QUERYNAME:   result.assign("queryname");  break;
        case HtsSortOrder::COORDINATE:  result.assign("coordinate"); break;
        default:                        result.assign("unknown");    break;
    }
    return result;
};
bool from_string(const char* value, HtsSortOrder& result) {
         if(value == NULL)                  result = HtsSortOrder::UNKNOWN;
    else if(!strcmp(value, "unsorted"))     result = HtsSortOrder::UNSORTED;
    else if(!strcmp(value, "queryname"))    result = HtsSortOrder::QUERYNAME;
    else if(!strcmp(value, "coordinate"))   result = HtsSortOrder::COORDINATE;
    else                                    result = HtsSortOrder::UNKNOWN;

    return (result == HtsSortOrder::UNKNOWN ? false : true);
};
bool from_string(const string& value, HtsSortOrder& result) {
    return from_string(value.c_str(), result);
};
void to_kstring(const HtsSortOrder& value, kstring_t& result) {
    ks_clear(result);
    string string_value(to_string(value));
    ks_put_string(string_value.c_str(), string_value.size(), result);
};
ostream& operator<<(ostream& o, const HtsSortOrder& value) {
    o << to_string(value);
    return o;
};
void encode_key_value(const string& key, const HtsSortOrder& value, Value& container, Document& document) {
    string string_value(to_string(value));
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

string to_string(const HtsGrouping& value) {
    string result;
    switch(value) {
        case HtsGrouping::QUERY:        result.assign("query");      break;
        case HtsGrouping::REFERENCE:    result.assign("reference");  break;
        default:                        result.assign("none");       break;
    }
    return result;
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
    string string_value(to_string(value));
    ks_put_string(string_value.c_str(), string_value.size(), result);
};
bool from_string(const string& value, HtsGrouping& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const HtsGrouping& value) {
    o << to_string(value);
    return o;
};
void encode_key_value(const string& key, const HtsGrouping& value, Value& container, Document& document) {
    string string_value(to_string(value));
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

string to_string(const htsFormatCategory& value) {
    string result;
    switch(value) {
        case htsFormatCategory::sequence_data:  result.assign("sequence data"); break;
        case htsFormatCategory::variant_data:   result.assign("variant data");  break;
        case htsFormatCategory::index_file:     result.assign("index file");    break;
        case htsFormatCategory::region_list:    result.assign("region list");   break;
        default:                                result.assign("unknown");       break;
    }
    return result;
};
bool from_string(const char* value, htsFormatCategory& result) {
         if(value == NULL)                      result = htsFormatCategory::unknown_category;
    else if(!strcmp(value, "sequence data"))    result = htsFormatCategory::sequence_data;
    else if(!strcmp(value, "variant data"))     result = htsFormatCategory::variant_data;
    else if(!strcmp(value, "index file"))       result = htsFormatCategory::index_file;
    else if(!strcmp(value, "region list"))      result = htsFormatCategory::region_list;
    else                                        result = htsFormatCategory::unknown_category;
    return (result == htsFormatCategory::unknown_category ? false : true);
};
ostream& operator<<(ostream& o, const htsFormatCategory& value) {
    o << to_string(value);
    return o;
};

string to_string(const htsExactFormat& value) {
    string result;
    switch(value) {
        case htsExactFormat::binary_format: result.assign("binary format"); break;
        case htsExactFormat::text_format:   result.assign("text format");   break;
        case htsExactFormat::sam:           result.assign("sam");           break;
        case htsExactFormat::bam:           result.assign("bam");           break;
        case htsExactFormat::bai:           result.assign("bai");           break;
        case htsExactFormat::cram:          result.assign("cram");          break;
        case htsExactFormat::crai:          result.assign("crai");          break;
        case htsExactFormat::vcf:           result.assign("vcf");           break;
        case htsExactFormat::bcf:           result.assign("bcf");           break;
        case htsExactFormat::csi:           result.assign("csi");           break;
        case htsExactFormat::gzi:           result.assign("gzi");           break;
        case htsExactFormat::tbi:           result.assign("tbi");           break;
        case htsExactFormat::bed:           result.assign("bed");           break;
        case htsExactFormat::htsget:        result.assign("htsget");        break;
        default:                            result.assign("unknown");       break;
    }
    return result;
};
bool from_string(const char* value, htsExactFormat& result) {
         if(value == NULL)                      result = htsExactFormat::unknown_format;
    else if(!strcmp(value, "binary format"))    result = htsExactFormat::binary_format;
    else if(!strcmp(value, "text format"))      result = htsExactFormat::text_format;
    else if(!strcmp(value, "sam"))              result = htsExactFormat::sam;
    else if(!strcmp(value, "bam"))              result = htsExactFormat::bam;
    else if(!strcmp(value, "bai"))              result = htsExactFormat::bai;
    else if(!strcmp(value, "cram"))             result = htsExactFormat::cram;
    else if(!strcmp(value, "crai"))             result = htsExactFormat::crai;
    else if(!strcmp(value, "vcf"))              result = htsExactFormat::vcf;
    else if(!strcmp(value, "bcf"))              result = htsExactFormat::bcf;
    else if(!strcmp(value, "csi"))              result = htsExactFormat::csi;
    else if(!strcmp(value, "gzi"))              result = htsExactFormat::gzi;
    else if(!strcmp(value, "tbi"))              result = htsExactFormat::tbi;
    else if(!strcmp(value, "bed"))              result = htsExactFormat::bed;
    else if(!strcmp(value, "htsget"))           result = htsExactFormat::htsget;
    else                                        result = htsExactFormat::unknown_format;
    return (result == htsExactFormat::unknown_format ? false : true);
};
ostream& operator<<(ostream& o, const htsExactFormat& value) {
    o << to_string(value);
    return o;
};

string to_string(const htsCompression& value) {
    string result;
    switch(value) {
        case htsCompression::gzip:      result.assign("gzip");              break;
        case htsCompression::bgzf:      result.assign("bgzf");              break;
        case htsCompression::custom:    result.assign("custom");            break;
        default:                        result.assign("no compression");    break;
    }
    return result;
};
bool from_string(const char* value, htsCompression& result) {
         if(value == NULL)              result = htsCompression::no_compression;
    else if(!strcmp(value, "gzip"))     result = htsCompression::gzip;
    else if(!strcmp(value, "bgzf"))     result = htsCompression::bgzf;
    else if(!strcmp(value, "custom"))   result = htsCompression::custom;
    else                                result = htsCompression::no_compression;
    return (result == htsCompression::no_compression ? false : true);
};
ostream& operator<<(ostream& o, const htsCompression& value) {
    o << to_string(value);
    return o;
};

/*
typedef struct htsFormat {
    enum htsFormatCategory category;
    enum htsExactFormat format;
    struct { short major, minor; } version;
    enum htsCompression compression;
    short compression_level;  // currently unused
    void *specific;  // format specific options; see struct hts_opt.
} htsFormat;
*/
ostream& operator<<(ostream& o, const htsFormat& hts_format) {
    o << "htsFormatCategory " << hts_format.category << endl;
    o << "htsExactFormat " << hts_format.format << endl;
    o << "version " << hts_format.version.major << "." << hts_format.version.minor << endl;
    o << "htsCompression " << hts_format.compression << endl;
    o << "compression_level " << int32_t(hts_format.compression_level) << endl;
    return o;
};


string to_string(const Platform& value) {
    string result;
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
    return result;
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
    string string_value(to_string(value));
    ks_put_string(string_value.c_str(), string_value.size(), result);
};
bool from_string(const string& value, Platform& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const Platform& value) {
    o << to_string(value);
    return o;
};
void encode_key_value(const string& key, const Platform& value, Value& container, Document& document) {
    string string_value(to_string(value));
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

string to_string(const Algorithm& value) {
    string result;
    switch(value) {
        case Algorithm::UNKNOWN:      result.assign("unknown");     break;
        case Algorithm::MDD:          result.assign("mdd");         break;
        case Algorithm::PAMLD:        result.assign("pamld");       break;
        case Algorithm::NAIVE:        result.assign("naive");       break;
        case Algorithm::PASSTHROUGH:  result.assign("passthrough"); break;
        case Algorithm::BENCHMARK:    result.assign("benchmark");   break;
        default:                                                    break;
    }
    return result;
};
bool from_string(const char* value, Algorithm& result) {
         if(value == NULL)                  result = Algorithm::UNKNOWN;
    else if(!strcmp(value, "mdd"))          result = Algorithm::MDD;
    else if(!strcmp(value, "pamld"))        result = Algorithm::PAMLD;
    else if(!strcmp(value, "naive"))        result = Algorithm::NAIVE;
    else if(!strcmp(value, "passthrough"))  result = Algorithm::PASSTHROUGH;
    else if(!strcmp(value, "benchmark"))    result = Algorithm::BENCHMARK;
    else                                    result = Algorithm::UNKNOWN;

    return (result == Algorithm::UNKNOWN ? false : true);
};
void to_kstring(const Algorithm& value, kstring_t& result) {
    ks_clear(result);
    string string_value(to_string(value));
    ks_put_string(string_value.c_str(), string_value.size(), result);
};
bool from_string(const string& value, Algorithm& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const Algorithm& value) {
    o << to_string(value);
    return o;
};
bool encode_key_value(const string& key, const Algorithm& value, Value& container, Document& document) {
    if(value != Algorithm::UNKNOWN) {
        string string_value(to_string(value));
        Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.RemoveMember(key.c_str());
        container.AddMember(k.Move(), v.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
template<> bool decode_value_by_key< Algorithm >(const Value::Ch* key, Algorithm& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
template <> Algorithm decode_value_by_key(const Value::Ch* key, const Value& container) {
    Algorithm value(Algorithm::UNKNOWN);
    decode_value_by_key(key, value, container);
    return value;
};

/* @HD The header line */
HeadHDAtom::HeadHDAtom() :
    VN({ 0, 0, NULL }),
    SO({ 0, 0, NULL }),
    GO({ 0, 0, NULL }),
    SS({ 0, 0, NULL }) {
    to_kstring(HtsSortOrder::UNKNOWN, SO);
    to_kstring(HtsGrouping::QUERY, GO);
};
HeadHDAtom::HeadHDAtom(const Value& ontology) try :
    VN(decode_value_by_key< kstring_t >("VN", ontology)),
    SO({ 0, 0, NULL }),
    GO({ 0, 0, NULL }),
    SS({ 0, 0, NULL }) {

    decode_value_by_key< kstring_t >("SO", SO, ontology);
    decode_value_by_key< kstring_t >("GO", GO, ontology);
    decode_value_by_key< kstring_t >("SS", SS, ontology);

    } catch(Error& error) {
        error.push("HeadHDAtom");
        throw;
};
HeadHDAtom::HeadHDAtom(const HeadHDAtom& other) :
    VN({ 0, 0, NULL }),
    SO({ 0, 0, NULL }),
    GO({ 0, 0, NULL }),
    SS({ 0, 0, NULL }) {

    if(ks_not_empty(other.VN)) ks_put_string(other.VN, VN);
    if(ks_not_empty(other.SO)) ks_put_string(other.SO, SO);
    if(ks_not_empty(other.GO)) ks_put_string(other.GO, GO);
    if(ks_not_empty(other.SS)) ks_put_string(other.SS, SS);
};
HeadHDAtom::~HeadHDAtom() {
    ks_free(VN);
    ks_free(SO);
    ks_free(GO);
    ks_free(SS);
};
HeadHDAtom& HeadHDAtom::operator=(const HeadHDAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(VN);
        ks_clear(SO);
        ks_clear(GO);
        ks_clear(SS);
        if(ks_not_empty(other.VN)) ks_put_string(other.VN, VN);
        if(ks_not_empty(other.SO)) ks_put_string(other.SO, SO);
        if(ks_not_empty(other.GO)) ks_put_string(other.GO, GO);
        if(ks_not_empty(other.SS)) ks_put_string(other.SS, SS);
    }
    return *this;
};
void HeadHDAtom::encode(kstring_t& buffer) const {
    ks_put_string_("@HD", 3, buffer);
    if(ks_not_empty(VN)) {
        ks_put_string_("\tVN:", 4, buffer);
        ks_put_string_(VN, buffer);
    }
    if(ks_not_empty(SO)) {
        ks_put_string_("\tSO:", 4, buffer);
        ks_put_string_(SO, buffer);
    }
    if(ks_not_empty(GO)) {
        ks_put_string_("\tGO:", 4, buffer);
        ks_put_string_(GO, buffer);
    }
    if(ks_not_empty(SS)) {
        ks_put_string_("\tSS:", 4, buffer);
        ks_put_string_(SS, buffer);
    }
    ks_put_character(LINE_BREAK, buffer);
};
char* HeadHDAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        ++position;
        uint16_t tag(tag_to_code(position));
        position += 3;
        switch (tag) {
            case uint16_t(HtsTagCode::VN): {
                position = copy_until_tag_end(position, end, VN);
                break;
            };
            case uint16_t(HtsTagCode::SO): {
                position = copy_until_tag_end(position, end, SO);
                break;
            };
            case uint16_t(HtsTagCode::GO): {
                position = copy_until_tag_end(position, end, GO);
                break;
            };
            case uint16_t(HtsTagCode::SS): {
                position = copy_until_tag_end(position, end, SS);
                break;
            };
            default:
                position = skip_to_tab(position, end);
                break;
        }
    }
    return ++position;
};
void HeadHDAtom::set_alignment_sort_order(const HtsSortOrder& order) {
    to_kstring(order, SO);
};
void HeadHDAtom::set_alignment_grouping(const HtsGrouping& grouping) {
    to_kstring(grouping, GO);
};
void HeadHDAtom::set_version(const htsFormat& format) {
    ks_clear(VN);
    ks_put_int32(format.version.major, VN);
    ks_put_character('.', VN);
    ks_put_int32(format.version.minor, VN);
};
ostream& operator<<(ostream& o, const HeadHDAtom& hd) {
    o << "@HD SAM head atom" << endl;
    if(ks_not_empty(hd.VN)) o << "VN : " << hd.VN.s << endl;
    if(ks_not_empty(hd.SO)) o << "SO : " << hd.SO.s << endl;
    if(ks_not_empty(hd.GO)) o << "GO : " << hd.GO.s << endl;
    if(ks_not_empty(hd.SS)) o << "SS : " << hd.SS.s << endl;
    return o;
};

/* @SQ reference sequence dictionary */
HeadSQAtom::HeadSQAtom() :
    index(-1),
    SN({ 0, 0, NULL }),
    LN(0),
    AH({ 0, 0, NULL }),
    AN({ 0, 0, NULL }),
    AS({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    M5({ 0, 0, NULL }),
    SP({ 0, 0, NULL }),
    TP({ 0, 0, NULL }),
    UR({ 0, 0, NULL }){
    ks_terminate(SN);
};
HeadSQAtom::HeadSQAtom(const Value& ontology) try :
    index(-1),
    SN(decode_value_by_key< kstring_t >("SN", ontology)),
    LN(decode_value_by_key< int32_t >("LN", ontology)),
    AH({ 0, 0, NULL }),
    AN({ 0, 0, NULL }),
    AS({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    M5({ 0, 0, NULL }),
    SP({ 0, 0, NULL }),
    TP({ 0, 0, NULL }),
    UR({ 0, 0, NULL }){

    decode_value_by_key< int32_t >("index", index, ontology);
    decode_value_by_key< kstring_t >("AH", AH, ontology);
    decode_value_by_key< kstring_t >("AN", AN, ontology);
    decode_value_by_key< kstring_t >("AS", AS, ontology);
    decode_value_by_key< kstring_t >("DS", DS, ontology);
    decode_value_by_key< kstring_t >("M5", M5, ontology);
    decode_value_by_key< kstring_t >("SP", SP, ontology);
    decode_value_by_key< kstring_t >("TP", TP, ontology);
    decode_value_by_key< kstring_t >("UR", UR, ontology);

    } catch(Error& error) {
        error.push("HeadSQAtom");
        throw;
};
HeadSQAtom::HeadSQAtom(const HeadSQAtom& other) :
    index(-1),
    SN({ 0, 0, NULL }),
    LN(0),
    AH({ 0, 0, NULL }),
    AN({ 0, 0, NULL }),
    AS({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    M5({ 0, 0, NULL }),
    SP({ 0, 0, NULL }),
    TP({ 0, 0, NULL }),
    UR({ 0, 0, NULL }) {

    ks_terminate(SN);
    if(other.index >= 0)       index = other.index;
    if(ks_not_empty(other.SN)) ks_put_string(other.SN, SN);
    if(other.LN > 0)           LN = other.LN;
    if(ks_not_empty(other.AH)) ks_put_string(other.AH, AH);
    if(ks_not_empty(other.AN)) ks_put_string(other.AN, AN);
    if(ks_not_empty(other.AS)) ks_put_string(other.AS, AS);
    if(ks_not_empty(other.DS)) ks_put_string(other.DS, DS);
    if(ks_not_empty(other.M5)) ks_put_string(other.M5, M5);
    if(ks_not_empty(other.SP)) ks_put_string(other.SP, SP);
    if(ks_not_empty(other.TP)) ks_put_string(other.TP, TP);
    if(ks_not_empty(other.UR)) ks_put_string(other.UR, UR);
};
HeadSQAtom::~HeadSQAtom() {
    ks_free(SN);
    ks_free(AH);
    ks_free(AN);
    ks_free(AS);
    ks_free(DS);
    ks_free(M5);
    ks_free(SP);
    ks_free(TP);
    ks_free(UR);
};
HeadSQAtom& HeadSQAtom::operator=(const HeadSQAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        index = -1;
        ks_clear(SN);
        LN = 0;
        ks_clear(AH);
        ks_clear(AN);
        ks_clear(AS);
        ks_clear(DS);
        ks_clear(M5);
        ks_clear(SP);
        ks_clear(TP);
        ks_clear(UR);
        if(other.index >= 0)       index = other.index;
        if(ks_not_empty(other.SN)) ks_put_string(other.SN, SN);
        if(other.LN > 0)           LN = other.LN;
        if(ks_not_empty(other.AH)) ks_put_string(other.AH, AH);
        if(ks_not_empty(other.AN)) ks_put_string(other.AN, AN);
        if(ks_not_empty(other.AS)) ks_put_string(other.AS, AS);
        if(ks_not_empty(other.DS)) ks_put_string(other.DS, DS);
        if(ks_not_empty(other.M5)) ks_put_string(other.M5, M5);
        if(ks_not_empty(other.SP)) ks_put_string(other.SP, SP);
        if(ks_not_empty(other.TP)) ks_put_string(other.TP, TP);
        if(ks_not_empty(other.UR)) ks_put_string(other.UR, UR);
    }
    return *this;
};
HeadSQAtom::operator string() const {
    return string(SN.s, SN.l);
};
void HeadSQAtom::encode(kstring_t& buffer) const {
    ks_put_string_("@SQ", 3, buffer);
    if(ks_not_empty(SN)) {
        ks_put_string_("\tSN:", 4, buffer);
        ks_put_string_(SN, buffer);
    }
    if(LN > 0) {
        ks_put_string_("\tLN:", 4, buffer);
        ks_put_int32(LN, buffer);
    }
    if(ks_not_empty(AH)) {
        ks_put_string_("\tAH:", 4, buffer);
        ks_put_string_(AH, buffer);
    }
    if(ks_not_empty(AN)) {
        ks_put_string_("\tAN:", 4, buffer);
        ks_put_string_(AN, buffer);
    }
    if(ks_not_empty(AS)) {
        ks_put_string_("\tAS:", 4, buffer);
        ks_put_string_(AS, buffer);
    }
    if(ks_not_empty(DS)) {
        ks_put_string_("\tDS:", 4, buffer);
        ks_put_string_(DS, buffer);
    }
    if(ks_not_empty(M5)) {
        ks_put_string_("\tM5:", 4, buffer);
        ks_put_string_(M5, buffer);
    }
    if(ks_not_empty(SP)) {
        ks_put_string_("\tSP:", 4, buffer);
        ks_put_string_(SP, buffer);
    }
    if(ks_not_empty(TP)) {
        ks_put_string_("\tTP:", 4, buffer);
        ks_put_string_(TP, buffer);
    }
    if(ks_not_empty(UR)) {
        ks_put_string_("\tUR:", 4, buffer);
        ks_put_string_(UR, buffer);
    }
    ks_put_character(LINE_BREAK, buffer);
};
char* HeadSQAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        ++position;
        uint16_t tag(tag_to_code(position));
        position += 3;
        switch (tag) {
            case uint16_t(HtsTagCode::SN): {
                position = copy_until_tag_end(position, end, SN);
                break;
            };
            case uint16_t(HtsTagCode::LN): {
                LN = static_cast< int32_t >(strtol(position, &position, 10));
                break;
            };
            case uint16_t(HtsTagCode::AH): {
                position = copy_until_tag_end(position, end, AH);
                break;
            };
            case uint16_t(HtsTagCode::AN): {
                position = copy_until_tag_end(position, end, AN);
                break;
            };
            case uint16_t(HtsTagCode::AS): {
                position = copy_until_tag_end(position, end, AS);
                break;
            };
            case uint16_t(HtsTagCode::DS): {
                position = copy_until_tag_end(position, end, DS);
                break;
            };
            case uint16_t(HtsTagCode::M5): {
                position = copy_until_tag_end(position, end, M5);
                break;
            };
            case uint16_t(HtsTagCode::SP): {
                position = copy_until_tag_end(position, end, SP);
                break;
            };
            case uint16_t(HtsTagCode::TP): {
                position = copy_until_tag_end(position, end, TP);
                break;
            };
            case uint16_t(HtsTagCode::UR): {
                position = copy_until_tag_end(position, end, UR);
                break;
            };
            default: {
                position = skip_to_tab(position, end);
                break;
            }
        }
    }
    return ++position;
};
void HeadSQAtom::update(const HeadSQAtom& other) {
    if(&other != this) {
        if(ks_not_empty(other.AH)) ks_assign_string(other.AH, AH);
        if(ks_not_empty(other.AN)) ks_assign_string(other.AN, AN);
        if(ks_not_empty(other.AS)) ks_assign_string(other.AS, AS);
        if(ks_not_empty(other.DS)) ks_assign_string(other.DS, DS);
        if(ks_not_empty(other.M5)) ks_assign_string(other.M5, M5);
        if(ks_not_empty(other.SP)) ks_assign_string(other.SP, SP);
        if(ks_not_empty(other.TP)) ks_assign_string(other.TP, TP);
        if(ks_not_empty(other.UR)) ks_assign_string(other.UR, UR);
    }
};
ostream& operator<<(ostream& o, const HeadSQAtom& sq) {
    o << "@SQ SAM head atom" << endl;
    if(sq.index   >= 0)     o << "index : " << sq.index << endl;
    if(ks_not_empty(sq.SN)) o << "SN : " << sq.SN.s << endl;
    if(sq.LN   > 0)         o << "LN : " << sq.LN   << endl;
    if(ks_not_empty(sq.AH)) o << "AH : " << sq.AH.s << endl;
    if(ks_not_empty(sq.AN)) o << "AN : " << sq.AN.s << endl;
    if(ks_not_empty(sq.AS)) o << "AS : " << sq.AS.s << endl;
    if(ks_not_empty(sq.DS)) o << "DS : " << sq.DS.s << endl;
    if(ks_not_empty(sq.M5)) o << "M5 : " << sq.M5.s << endl;
    if(ks_not_empty(sq.SP)) o << "SP : " << sq.SP.s << endl;
    if(ks_not_empty(sq.TP)) o << "TP : " << sq.TP.s << endl;
    if(ks_not_empty(sq.UR)) o << "UR : " << sq.UR.s << endl;
    return o;
};

/* @RG Read Group */
HeadRGAtom::HeadRGAtom() :
    index(-1),
    ID({ 0, 0, NULL }),
    BC({ 0, 0, NULL }),
    CN({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    DT({ 0, 0, NULL }),
    FO({ 0, 0, NULL }),
    KS({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    PI({ 0, 0, NULL }),
    PL({ 0, 0, NULL }),
    PM({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    SM({ 0, 0, NULL }) {
    ks_terminate(ID);
};
HeadRGAtom::HeadRGAtom(const Value& ontology) try :
    index(-1),
    ID(decode_value_by_key< kstring_t >("ID", ontology)),
    BC({ 0, 0, NULL }),
    CN({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    DT({ 0, 0, NULL }),
    FO({ 0, 0, NULL }),
    KS({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    PI({ 0, 0, NULL }),
    PL({ 0, 0, NULL }),
    PM({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    SM({ 0, 0, NULL }) {

    decode_value_by_key< int32_t >("index", index, ontology);
    decode_value_by_key< kstring_t >("BC", BC, ontology);
    decode_value_by_key< kstring_t >("CN", CN, ontology);
    decode_value_by_key< kstring_t >("DS", DS, ontology);
    decode_value_by_key< kstring_t >("DT", DT, ontology);
    decode_value_by_key< kstring_t >("FO", FO, ontology);
    decode_value_by_key< kstring_t >("KS", KS, ontology);
    decode_value_by_key< kstring_t >("LB", LB, ontology);
    decode_value_by_key< kstring_t >("PG", PG, ontology);
    decode_value_by_key< kstring_t >("PI", PI, ontology);
    decode_value_by_key< kstring_t >("PL", PL, ontology);
    decode_value_by_key< kstring_t >("PM", PM, ontology);
    decode_value_by_key< kstring_t >("PU", PU, ontology);
    decode_value_by_key< kstring_t >("SM", SM, ontology);

    } catch(Error& error) {
        error.push("HeadRGAtom");
        throw;
};
HeadRGAtom::HeadRGAtom(const HeadRGAtom& other) :
    index(-1),
    ID({ 0, 0, NULL }),
    BC({ 0, 0, NULL }),
    CN({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    DT({ 0, 0, NULL }),
    FO({ 0, 0, NULL }),
    KS({ 0, 0, NULL }),
    LB({ 0, 0, NULL }),
    PG({ 0, 0, NULL }),
    PI({ 0, 0, NULL }),
    PL({ 0, 0, NULL }),
    PM({ 0, 0, NULL }),
    PU({ 0, 0, NULL }),
    SM({ 0, 0, NULL }) {

    ks_terminate(ID);
    if(other.index >= 0)       index = other.index;
    if(ks_not_empty(other.ID)) ks_put_string(other.ID, ID);
    if(ks_not_empty(other.BC)) ks_put_string(other.BC, BC);
    if(ks_not_empty(other.CN)) ks_put_string(other.CN, CN);
    if(ks_not_empty(other.DS)) ks_put_string(other.DS, DS);
    if(ks_not_empty(other.DT)) ks_put_string(other.DT, DT);
    if(ks_not_empty(other.FO)) ks_put_string(other.FO, FO);
    if(ks_not_empty(other.KS)) ks_put_string(other.KS, KS);
    if(ks_not_empty(other.LB)) ks_put_string(other.LB, LB);
    if(ks_not_empty(other.PG)) ks_put_string(other.PG, PG);
    if(ks_not_empty(other.PI)) ks_put_string(other.PI, PI);
    if(ks_not_empty(other.PL)) ks_put_string(other.PL, PL);
    if(ks_not_empty(other.PM)) ks_put_string(other.PM, PM);
    if(ks_not_empty(other.PU)) ks_put_string(other.PU, PU);
    if(ks_not_empty(other.SM)) ks_put_string(other.SM, SM);
};
HeadRGAtom::~HeadRGAtom() {
    ks_free(ID);
    ks_free(BC);
    ks_free(CN);
    ks_free(DS);
    ks_free(DT);
    ks_free(FO);
    ks_free(KS);
    ks_free(LB);
    ks_free(PG);
    ks_free(PI);
    ks_free(PL);
    ks_free(PM);
    ks_free(PU);
    ks_free(SM);
};
HeadRGAtom& HeadRGAtom::operator=(const HeadRGAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        index = -1;
        ks_clear(ID);
        ks_clear(BC);
        ks_clear(CN);
        ks_clear(DS);
        ks_clear(DT);
        ks_clear(FO);
        ks_clear(KS);
        ks_clear(LB);
        ks_clear(PG);
        ks_clear(PI);
        ks_clear(PL);
        ks_clear(PM);
        ks_clear(PU);
        ks_clear(SM);
        if(other.index >= 0)       index = other.index;
        if(ks_not_empty(other.ID)) ks_put_string(other.ID, ID);
        if(ks_not_empty(other.BC)) ks_put_string(other.BC, BC);
        if(ks_not_empty(other.CN)) ks_put_string(other.CN, CN);
        if(ks_not_empty(other.DS)) ks_put_string(other.DS, DS);
        if(ks_not_empty(other.DT)) ks_put_string(other.DT, DT);
        if(ks_not_empty(other.FO)) ks_put_string(other.FO, FO);
        if(ks_not_empty(other.KS)) ks_put_string(other.KS, KS);
        if(ks_not_empty(other.LB)) ks_put_string(other.LB, LB);
        if(ks_not_empty(other.PG)) ks_put_string(other.PG, PG);
        if(ks_not_empty(other.PI)) ks_put_string(other.PI, PI);
        if(ks_not_empty(other.PL)) ks_put_string(other.PL, PL);
        if(ks_not_empty(other.PM)) ks_put_string(other.PM, PM);
        if(ks_not_empty(other.PU)) ks_put_string(other.PU, PU);
        if(ks_not_empty(other.SM)) ks_put_string(other.SM, SM);
    }
    return *this;
};
HeadRGAtom::operator string() const {
    return string(ID.s, ID.l);
};
void HeadRGAtom::encode(kstring_t& buffer) const {
    ks_put_string_("@RG", 3, buffer);
    if(ks_not_empty(ID)) {
        ks_put_string_("\tID:", 4, buffer);
        ks_put_string_(ID, buffer);
    }
    if(ks_not_empty(BC)) {
        ks_put_string_("\tBC:", 4, buffer);
        ks_put_string_(BC, buffer);
    }
    if(ks_not_empty(CN)) {
        ks_put_string_("\tCN:", 4, buffer);
        ks_put_string_(CN, buffer);
    }
    if(ks_not_empty(DS)) {
        ks_put_string_("\tDS:", 4, buffer);
        ks_put_string_(DS, buffer);
    }
    if(ks_not_empty(DT)) {
        ks_put_string_("\tDT:", 4, buffer);
        ks_put_string_(DT, buffer);
    }
    if(ks_not_empty(FO)) {
        ks_put_string_("\tFO:", 4, buffer);
        ks_put_string_(FO, buffer);
    }
    if(ks_not_empty(KS)) {
        ks_put_string_("\tKS:", 4, buffer);
        ks_put_string_(KS, buffer);
    }
    if(ks_not_empty(LB)) {
        ks_put_string_("\tLB:", 4, buffer);
        ks_put_string_(LB, buffer);
    }
    if(ks_not_empty(PG)) {
        ks_put_string_("\tPG:", 4, buffer);
        ks_put_string_(PG, buffer);
    }
    if(ks_not_empty(PI)) {
        ks_put_string_("\tPI:", 4, buffer);
        ks_put_string_(PI, buffer);
    }
    if(ks_not_empty(PL)) {
        ks_put_string_("\tPL:", 4, buffer);
        ks_put_string_(PL, buffer);
    }
    if(ks_not_empty(PM)) {
        ks_put_string_("\tPM:", 4, buffer);
        ks_put_string_(PM, buffer);
    }
    if(ks_not_empty(PU)) {
        ks_put_string_("\tPU:", 4, buffer);
        ks_put_string_(PU, buffer);
    }
    if(ks_not_empty(SM)) {
        ks_put_string_("\tSM:", 4, buffer);
        ks_put_string_(SM, buffer);
    }
    ks_put_character(LINE_BREAK, buffer);
};
char* HeadRGAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        ++position;
        uint16_t tag(tag_to_code(position));
        position += 3;
        switch (tag) {
            case uint16_t(HtsTagCode::ID): {
                position = copy_until_tag_end(position, end, ID);
                break;
            };
            case uint16_t(HtsTagCode::BC): {
                position = copy_until_tag_end(position, end, BC);
                break;
            };
            case uint16_t(HtsTagCode::CN): {
                position = copy_until_tag_end(position, end, CN);
                break;
            };
            case uint16_t(HtsTagCode::DS): {
                position = copy_until_tag_end(position, end, DS);
                break;
            };
            case uint16_t(HtsTagCode::DT): {
                position = copy_until_tag_end(position, end, DT);
                break;
            };
            case uint16_t(HtsTagCode::FO): {
                position = copy_until_tag_end(position, end, FO);
                break;
            };
            case uint16_t(HtsTagCode::KS): {
                position = copy_until_tag_end(position, end, KS);
                break;
            };
            case uint16_t(HtsTagCode::LB): {
                position = copy_until_tag_end(position, end, LB);
                break;
            };
            case uint16_t(HtsTagCode::PG): {
                position = copy_until_tag_end(position, end, PG);
                break;
            };
            case uint16_t(HtsTagCode::PI): {
                position = copy_until_tag_end(position, end, PI);
                break;
            };
            case uint16_t(HtsTagCode::PL): {
                position = copy_until_tag_end(position, end, PL);
                break;
            };
            case uint16_t(HtsTagCode::PM): {
                position = copy_until_tag_end(position, end, PM);
                break;
            };
            case uint16_t(HtsTagCode::PU): {
                position = copy_until_tag_end(position, end, PU);
                break;
            };
            case uint16_t(HtsTagCode::SM): {
                position = copy_until_tag_end(position, end, SM);
                break;
            };
            default:
                position = skip_to_tab(position, end);
                break;
        }
    }
    return ++position;
};
void HeadRGAtom::update(const HeadRGAtom& other) {
    if(&other != this) {
        if(ks_not_empty(other.BC)) ks_assign_string(other.BC, BC);
        if(ks_not_empty(other.CN)) ks_assign_string(other.CN, CN);
        if(ks_not_empty(other.DS)) ks_assign_string(other.DS, DS);
        if(ks_not_empty(other.DT)) ks_assign_string(other.DT, DT);
        if(ks_not_empty(other.FO)) ks_assign_string(other.FO, FO);
        if(ks_not_empty(other.KS)) ks_assign_string(other.KS, KS);
        if(ks_not_empty(other.LB)) ks_assign_string(other.LB, LB);
        if(ks_not_empty(other.PG)) ks_assign_string(other.PG, PG);
        if(ks_not_empty(other.PI)) ks_assign_string(other.PI, PI);
        if(ks_not_empty(other.PL)) ks_assign_string(other.PL, PL);
        if(ks_not_empty(other.PM)) ks_assign_string(other.PM, PM);
        if(ks_not_empty(other.PU)) ks_assign_string(other.PU, PU);
        if(ks_not_empty(other.SM)) ks_assign_string(other.SM, SM);
    }
};
void HeadRGAtom::set_platform(const Platform& value) {
    to_kstring(value, PL);
};
ostream& operator<<(ostream& o, const HeadRGAtom& rg) {
    o << "@RG SAM head atom" << endl;
    if(rg.index   >= 0)     o << "index : " << rg.index << endl;
    if(ks_not_empty(rg.ID)) o << "ID : " << rg.ID.s << endl;
    if(ks_not_empty(rg.BC)) o << "BC : " << rg.BC.s << endl;
    if(ks_not_empty(rg.CN)) o << "CN : " << rg.CN.s << endl;
    if(ks_not_empty(rg.DS)) o << "DS : " << rg.DS.s << endl;
    if(ks_not_empty(rg.DT)) o << "DT : " << rg.DT.s << endl;
    if(ks_not_empty(rg.FO)) o << "FO : " << rg.FO.s << endl;
    if(ks_not_empty(rg.KS)) o << "KS : " << rg.KS.s << endl;
    if(ks_not_empty(rg.LB)) o << "LB : " << rg.LB.s << endl;
    if(ks_not_empty(rg.PG)) o << "PG : " << rg.PG.s << endl;
    if(ks_not_empty(rg.PI)) o << "PI : " << rg.PI.s << endl;
    if(ks_not_empty(rg.PL)) o << "PL : " << rg.PL.s << endl;
    if(ks_not_empty(rg.PM)) o << "PM : " << rg.PM.s << endl;
    if(ks_not_empty(rg.PU)) o << "PU : " << rg.PU.s << endl;
    if(ks_not_empty(rg.SM)) o << "SM : " << rg.SM.s << endl;
    return o;
};
template<> bool decode_value_by_key< vector< HeadRGAtom > >(const Value::Ch* key, vector< HeadRGAtom >& value, const Value& container) {
    bool result(false);
    if(container.IsObject()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd() && !reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                value.clear();
                Value::ConstMemberIterator undetermined_reference = reference->value.FindMember("undetermined");
                if(undetermined_reference != reference->value.MemberEnd()) {
                    result = true;
                    Value::ConstMemberIterator codec_reference = reference->value.FindMember("codec");
                    if(codec_reference != reference->value.MemberEnd()) {
                        value.reserve(codec_reference->value.MemberCount() + 1);
                        HeadRGAtom undetermined_rg(undetermined_reference->value);
                        value.emplace_back(undetermined_rg);
                        for(auto& record : codec_reference->value.GetObject()) {
                            HeadRGAtom rg(record.value);
                            value.emplace_back(rg);
                        }
                    } else {
                        value.reserve(1);
                        HeadRGAtom undetermined_rg(undetermined_reference->value);
                        value.emplace_back(undetermined_rg);
                    }
                } else { throw ConfigurationError("classifier must declare an undetermined element"); }
            }
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return result;
};
bool encode_value(const HeadRGAtom& value, Value& container, Document& document) {
    if(container.IsObject()) {
        encode_key_value("ID", value.ID, container, document);
        encode_key_value("BC", value.BC, container, document);
        encode_key_value("CN", value.CN, container, document);
        encode_key_value("DS", value.DS, container, document);
        encode_key_value("DT", value.DT, container, document);
        encode_key_value("FO", value.FO, container, document);
        encode_key_value("KS", value.KS, container, document);
        encode_key_value("LB", value.LB, container, document);
        encode_key_value("PG", value.PG, container, document);
        encode_key_value("PI", value.PI, container, document);
        encode_key_value("PL", value.PL, container, document);
        encode_key_value("PM", value.PM, container, document);
        encode_key_value("PU", value.PU, container, document);
        encode_key_value("SM", value.SM, container, document);
        return true;
    } else { throw ConfigurationError("Read Group element must be a dictionary"); }
};
bool encode_key_value(const string& key, const list< HeadRGAtom >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            Value c(kObjectType);
            encode_value(v, c, document);
            array.PushBack(c.Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};

/* @PG program */
HeadPGAtom::HeadPGAtom() :
    index(-1),
    ID({ 0, 0, NULL }),
    PN({ 0, 0, NULL }),
    CL({ 0, 0, NULL }),
    PP({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    VN({ 0, 0, NULL }){
    ks_terminate(ID);
};
HeadPGAtom::HeadPGAtom(const Value& ontology) try :
    index(-1),
    ID(decode_value_by_key< kstring_t >("ID", ontology)),
    PN({ 0, 0, NULL }),
    CL({ 0, 0, NULL }),
    PP({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    VN({ 0, 0, NULL }){

    decode_value_by_key< int32_t >("index", index, ontology);
    decode_value_by_key< kstring_t >("PN", PN, ontology);
    decode_value_by_key< kstring_t >("CL", CL, ontology);
    decode_value_by_key< kstring_t >("PP", PP, ontology);
    decode_value_by_key< kstring_t >("DS", DS, ontology);
    decode_value_by_key< kstring_t >("VN", VN, ontology);

    } catch(Error& error) {
        error.push("HeadPGAtom");
        throw;
};
HeadPGAtom::HeadPGAtom(const HeadPGAtom& other) :
    index(-1),
    ID({ 0, 0, NULL }),
    PN({ 0, 0, NULL }),
    CL({ 0, 0, NULL }),
    PP({ 0, 0, NULL }),
    DS({ 0, 0, NULL }),
    VN({ 0, 0, NULL }){

    ks_terminate(ID);
    if(other.index >= 0)       index = other.index;
    if(ks_not_empty(other.ID)) ks_put_string(other.ID, ID);
    if(ks_not_empty(other.PN)) ks_put_string(other.PN, PN);
    if(ks_not_empty(other.CL)) ks_put_string(other.CL, CL);
    if(ks_not_empty(other.PP)) ks_put_string(other.PP, PP);
    if(ks_not_empty(other.DS)) ks_put_string(other.DS, DS);
    if(ks_not_empty(other.VN)) ks_put_string(other.VN, VN);
};
HeadPGAtom::~HeadPGAtom() {
    ks_free(ID);
    ks_free(PN);
    ks_free(CL);
    ks_free(PP);
    ks_free(DS);
    ks_free(VN);
};
HeadPGAtom& HeadPGAtom::operator=(const HeadPGAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        index = -1;
        ks_clear(ID);
        ks_clear(PN);
        ks_clear(CL);
        ks_clear(PP);
        ks_clear(DS);
        ks_clear(VN);
        if(other.index >= 0)       index = other.index;
        if(ks_not_empty(other.ID)) ks_put_string(other.ID, ID);
        if(ks_not_empty(other.PN)) ks_put_string(other.PN, PN);
        if(ks_not_empty(other.CL)) ks_put_string(other.CL, CL);
        if(ks_not_empty(other.PP)) ks_put_string(other.PP, PP);
        if(ks_not_empty(other.DS)) ks_put_string(other.DS, DS);
        if(ks_not_empty(other.VN)) ks_put_string(other.VN, VN);
    }
    return *this;
};
HeadPGAtom::operator string() const {
    return string(ID.s, ID.l);
};
void HeadPGAtom::encode(kstring_t& buffer) const {
    ks_put_string_("@PG", 3, buffer);
    if(ks_not_empty(ID)) {
        ks_put_string_("\tID:", 4, buffer);
        ks_put_string_(ID, buffer);
    }
    if(ks_not_empty(PN)) {
        ks_put_string_("\tPN:", 4, buffer);
        ks_put_string_(PN, buffer);
    }
    if(ks_not_empty(CL)) {
        ks_put_string_("\tCL:", 4, buffer);
        ks_put_string_(CL, buffer);
    }
    if(ks_not_empty(PP)) {
        ks_put_string_("\tPP:", 4, buffer);
        ks_put_string_(PP, buffer);
    }
    if(ks_not_empty(DS)) {
        ks_put_string_("\tDS:", 4, buffer);
        ks_put_string_(DS, buffer);
    }
    if(ks_not_empty(VN)) {
        ks_put_string_("\tVN:", 4, buffer);
        ks_put_string_(VN, buffer);
    }
    ks_put_character(LINE_BREAK, buffer);
};
char* HeadPGAtom::decode(char* position, const char* end) {
    while(*position == '\t' && position <= end) {
        ++position;
        uint16_t tag(tag_to_code(position));
        position += 3;
        switch (tag) {
            case uint16_t(HtsTagCode::ID): {
                position = copy_until_tag_end(position, end, ID);
                break;
            };
            case uint16_t(HtsTagCode::PN): {
                position = copy_until_tag_end(position, end, PN);
                break;
            };
            case uint16_t(HtsTagCode::CL): {
                position = copy_until_tag_end(position, end, CL);
                break;
            };
            case uint16_t(HtsTagCode::PP): {
                position = copy_until_tag_end(position, end, PP);
                break;
            };
            case uint16_t(HtsTagCode::DS): {
                position = copy_until_tag_end(position, end, DS);
                break;
            };
            case uint16_t(HtsTagCode::VN): {
                position = copy_until_tag_end(position, end, VN);
                break;
            };
            default: {
                position = skip_to_tab(position, end);
                break;
            }
        }
    }
    return ++position;
};
void HeadPGAtom::update(const HeadPGAtom& other) {
    if(&other != this) {
        if(ks_not_empty(other.PN)) ks_assign_string(other.PN, PN);
        if(ks_not_empty(other.CL)) ks_assign_string(other.CL, CL);
        if(ks_not_empty(other.PP)) ks_assign_string(other.PP, PP);
        if(ks_not_empty(other.DS)) ks_assign_string(other.DS, DS);
        if(ks_not_empty(other.VN)) ks_assign_string(other.VN, VN);
    }
};
ostream& operator<<(ostream& o, const HeadPGAtom& pg) {
    o << "@PG SAM head atom" << endl;
    if(pg.index   >= 0)     o << "index : " << pg.index << endl;
    if(ks_not_empty(pg.ID)) o << "ID : " << pg.ID.s << endl;
    if(ks_not_empty(pg.PN)) o << "PN : " << pg.PN.s << endl;
    if(ks_not_empty(pg.CL)) o << "CL : " << pg.CL.s << endl;
    if(ks_not_empty(pg.PP)) o << "PP : " << pg.PP.s << endl;
    if(ks_not_empty(pg.DS)) o << "DS : " << pg.DS.s << endl;
    if(ks_not_empty(pg.VN)) o << "VN : " << pg.VN.s << endl;
    return o;
};
template<> HeadPGAtom decode_value_by_key< HeadPGAtom >(const Value::Ch* key, const Value& container) {
    if(container.IsObject()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            HeadPGAtom PG(reference->value);
            return PG;
        } else { throw ConfigurationError(string(key) + " not found"); }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
};

/* @CO free text comment */
HeadCOAtom::HeadCOAtom() :
    CO({ 0, 0, NULL }) {
};
HeadCOAtom::HeadCOAtom(const Value& ontology) try :
    CO({ 0, 0, NULL }) {

    decode_value_by_key< kstring_t >("CO", CO, ontology);

    } catch(Error& error) {
        error.push("HeadCOAtom");
        throw;
};
HeadCOAtom::HeadCOAtom(const HeadCOAtom& other) :
    CO({ 0, 0, NULL }){
    if(ks_not_empty(other.CO)) ks_put_string(other.CO, CO);
};
HeadCOAtom::~HeadCOAtom() {
    ks_free(CO);
};
HeadCOAtom& HeadCOAtom::operator=(const HeadCOAtom& other) {
    if(&other == this) {
        return *this;
    } else {
        ks_clear(CO);
        if(ks_not_empty(other.CO)) ks_put_string(other.CO, CO);
    }
    return *this;
};
char* HeadCOAtom::decode(char* position, const char* end) {
    if(*position == '\t' && position <= end) {
        ++position;
        position = copy_until_linebreak(position, end, CO);
    }
    return ++position;
};
void HeadCOAtom::encode(kstring_t& buffer) const {
    if(ks_not_empty(CO)) {
        ks_put_string_("@CO:", 4, buffer);
        ks_put_string_(CO, buffer);
        ks_put_character(LINE_BREAK, buffer);
    }
};
ostream& operator<<(ostream& o, const HeadCOAtom& co) {
    o << "@CO SAM head atom" << endl;
    if(ks_not_empty(co.CO)) o << "CO : " << co.CO.s << endl;
    return o;
};

static inline char* skip_to_linebreak(char* source, const char* end) {
    char* position(source);
    while (*position != LINE_BREAK && position < end) {
        ++position;
    }
    return position;
};

void HtsHead::decode(const sam_hdr_t* hdr) {
    if(hdr != NULL) {
        /*  read in @SQ element
            Seems like if only the SN and LN atoms are present
            @SQ records are absent from plain text section */
        for(int32_t index(0); index < hdr->n_targets; ++index) {
            HeadSQAtom sq;
            sq.index = index;
            sq.LN = hdr->target_len[index];
            ks_put_string(hdr->target_name[index], sq.SN);
            add_reference(sq);
        }

        char* position(hdr->text);
        char* end(position + hdr->l_text);
        while(position < end) {
            if(*position == '@') {
                ++position;
                uint16_t code = tag_to_code(position);
                position += 2;
                switch (code) {
                    case uint16_t(HtsTagCode::HD): {
                        position = hd.decode(position, end);
                        break;
                    };
                    case uint16_t(HtsTagCode::RG): {
                        HeadRGAtom rg;
                        position = rg.decode(position, end);
                        add_read_group(rg);
                        break;
                    };
                    case uint16_t(HtsTagCode::PG): {
                        HeadPGAtom pg;
                        position = pg.decode(position, end);
                        add_program(pg);
                        break;
                    };
                    case uint16_t(HtsTagCode::CO): {
                        HeadCOAtom co;
                        position = co.decode(position, end);
                        add_comment(co);
                        break;
                    };
                    case uint16_t(HtsTagCode::SQ): {
                        HeadSQAtom sq;
                        position = sq.decode(position, end);
                        add_reference(sq);
                        break;
                    };
                    default:
                        position = skip_to_linebreak(position, end);
                        ++position;
                        break;
                }
            }
        }
    }
};
void HtsHead::encode(sam_hdr_t* hdr) const {
    kstring_t buffer = { 0, 0, NULL };
    hd.encode(buffer);

    vector< const HeadSQAtom* > reference_by_index(SQ_by_id.size());
    for(const auto& record : SQ_by_id) {
        reference_by_index[record.second.index] = &record.second;
    }
    for(const auto& record : reference_by_index) {
        record->encode(buffer);
    }

    vector< const HeadRGAtom* > read_group_by_index(RG_by_id.size());
    for(const auto& record : RG_by_id) {
        read_group_by_index[record.second.index] = &record.second;
    }
    for(const auto& record : read_group_by_index) {
        record->encode(buffer);
    }

    vector< const HeadPGAtom* > program_by_index(PG_by_id.size());
    for(const auto& record : PG_by_id) {
        program_by_index[record.second.index] = &record.second;
    }
    for(const auto& record : program_by_index) {
        record->encode(buffer);
    }

    for(const auto& record: comment_by_index){
        record.encode(buffer);
    }

    if(buffer.l <= numeric_limits< uint32_t >::max()) {
        hdr->n_targets = 0;
        hdr->l_text = static_cast< uint32_t >(buffer.l);
        if((hdr->text = static_cast< char* >(malloc(hdr->l_text + 1))) == NULL) {
            throw OutOfMemoryError();
        }
        memcpy(hdr->text, buffer.s, hdr->l_text + 1);
        ks_free(buffer);
    } else { throw OverflowError("SAM header must not exceed " + to_string(numeric_limits< uint32_t >::max()) + " bytes"); }
};
void HtsHead::add_reference(const HeadSQAtom& atom) {
    string key(atom);
    auto record = SQ_by_id.find(key);
    if(record != SQ_by_id.end()) {
        record->second.update(atom);
    } else {
        HeadSQAtom o(atom);
        o.index = SQ_by_id.size();
        SQ_by_id.insert(make_pair(key, o));
    }
};
void HtsHead::add_read_group(const HeadRGAtom& atom) {
    string key(atom);
    auto record = RG_by_id.find(key);
    if(record != RG_by_id.end()) {
        record->second.update(atom);
    } else {
        HeadRGAtom o(atom);
        o.index = RG_by_id.size();
        RG_by_id.insert(make_pair(key, o));
    }
};
void HtsHead::add_program(const HeadPGAtom& atom) {
    string key(atom);
    auto record = PG_by_id.find(key);
    if(record != PG_by_id.end()) {
        record->second.update(atom);
    } else {
        HeadPGAtom o(atom);
        o.index = PG_by_id.size();
        PG_by_id.insert(make_pair(key, o));
    }
};
void HtsHead::add_comment(const HeadCOAtom& atom) {
    comment_by_index.push_back(atom);
};
ostream& operator<<(ostream& o, const HtsHead& head) {
    o << head.hd;
    for(const auto& record : head.SQ_by_id) {
        o << record.second;
    }
    for(const auto& record : head.RG_by_id) {
        o << record.second;
    }
    for(const auto& record : head.PG_by_id) {
        o << record.second;
    }
    for(const auto& co : head.comment_by_index) {
        o << co;
    }
    return o;
};
