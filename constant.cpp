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

ostream& operator<<(ostream& o, const ProgramAction& type) {
    string string_value;
    string_value << type;
    o << string_value;
    return o;
};
string& operator<<(string& o, const ProgramAction& type) {
    switch (type) {
        case ProgramAction::DEMULTIPLEX:    o.assign("demux");      break;
        case ProgramAction::QUALITY:        o.assign("quality");    break;
        default:                            o.assign("unknown");    break;
    }
    return o;
};
void operator>>(const string& s, ProgramAction& type) {
    if(s == "demux")            type = ProgramAction::DEMULTIPLEX;
    else if(s == "quality")     type = ProgramAction::QUALITY;
    else                        type = ProgramAction::UNKNOWN;
};
void encode_key_value(const string& key, const ProgramAction& value, Value& container, Document& document) {
    string string_value;
    string_value << value;
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
void decode_program_action_by_key(const Value::Ch* key, ProgramAction& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsString()) {
            string string_value(element->value.GetString(), element->value.GetStringLength());
            string_value >> value;
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
};

ostream& operator<<(ostream& o, const FormatType& type) {
    switch (type) {
        case FormatType::FASTQ: o << "fastq";   break;
        case FormatType::SAM:   o << "sam";     break;
        case FormatType::BAM:   o << "bam";     break;
        case FormatType::BAI:   o << "bai";     break;
        case FormatType::CRAM:  o << "cram";    break;
        case FormatType::CRAI:  o << "crai";    break;
        case FormatType::VCF:   o << "vcf";     break;
        case FormatType::BCF:   o << "bcf";     break;
        case FormatType::CSI:   o << "csi";     break;
        case FormatType::GZI:   o << "gzi";     break;
        case FormatType::TBI:   o << "tbi";     break;
        case FormatType::BED:   o << "bed";     break;
        case FormatType::JSON:  o << "json";    break;
        default:                o << "unknown"; break;
    }
    return o;
};
string& operator<<(string& o, const FormatType& type) {
    switch (type) {
        case FormatType::FASTQ: o.assign("fastq");  break;
        case FormatType::SAM:   o.assign("sam");    break;
        case FormatType::BAM:   o.assign("bam");    break;
        case FormatType::BAI:   o.assign("bai");    break;
        case FormatType::CRAM:  o.assign("cram");   break;
        case FormatType::CRAI:  o.assign("crai");   break;
        case FormatType::VCF:   o.assign("vcf");    break;
        case FormatType::BCF:   o.assign("bcf");    break;
        case FormatType::CSI:   o.assign("csi");    break;
        case FormatType::GZI:   o.assign("gzi");    break;
        case FormatType::TBI:   o.assign("tbi");    break;
        case FormatType::BED:   o.assign("bed");    break;
        case FormatType::JSON:  o.assign("json");   break;
        default:                                    break;
    }
    return o;
};
void operator>>(const char* s, FormatType& type) {
         if(s == NULL)              type = FormatType::UNKNOWN;
    else if(!strcmp(s, "fastq"))    type = FormatType::FASTQ;
    else if(!strcmp(s, "sam"))      type = FormatType::SAM;
    else if(!strcmp(s, "bam"))      type = FormatType::BAM;
    else if(!strcmp(s, "bai"))      type = FormatType::BAI;
    else if(!strcmp(s, "cram"))     type = FormatType::CRAM;
    else if(!strcmp(s, "crai"))     type = FormatType::CRAI;
    else if(!strcmp(s, "vcf"))      type = FormatType::VCF;
    else if(!strcmp(s, "bcf"))      type = FormatType::BCF;
    else if(!strcmp(s, "csi"))      type = FormatType::CSI;
    else if(!strcmp(s, "gzi"))      type = FormatType::GZI;
    else if(!strcmp(s, "TBI"))      type = FormatType::TBI;
    else if(!strcmp(s, "bed"))      type = FormatType::BED;
    else if(!strcmp(s, "json"))     type = FormatType::JSON;
    else                            type = FormatType::UNKNOWN;
};
void operator>>(const string& s, FormatType& type) {
         if(s == "fastq")   type = FormatType::FASTQ;
    else if(s == "sam")     type = FormatType::SAM;
    else if(s == "bam")     type = FormatType::BAM;
    else if(s == "bai")     type = FormatType::BAI;
    else if(s == "cram")    type = FormatType::CRAM;
    else if(s == "crai")    type = FormatType::CRAI;
    else if(s == "vcf")     type = FormatType::VCF;
    else if(s == "bcf")     type = FormatType::BCF;
    else if(s == "csi")     type = FormatType::CSI;
    else if(s == "gzi")     type = FormatType::GZI;
    else if(s == "TBI")     type = FormatType::TBI;
    else if(s == "bed")     type = FormatType::BED;
    else if(s == "json")    type = FormatType::JSON;
    else                    type = FormatType::UNKNOWN;
}

ostream& operator<<(ostream& o, const HtsSortOrder& order) {
    switch (order) {
        case HtsSortOrder::UNSORTED:    o << "unsorted";    break;
        case HtsSortOrder::QUERYNAME:   o << "queryname";   break;
        case HtsSortOrder::COORDINATE:  o << "coordinate";  break;
        default:                        o << "unknown";     break;
    }
    return o;
};
string& operator<<(string& o, const HtsSortOrder& order) {
    switch (order) {
        case HtsSortOrder::UNSORTED:    o.assign("unsorted");   break;
        case HtsSortOrder::QUERYNAME:   o.assign("queryname");  break;
        case HtsSortOrder::COORDINATE:  o.assign("coordinate"); break;
        default:                        o.assign("unknown");    break;
    }
    return o;
};
kstring_t& operator<<(kstring_t& o, const HtsSortOrder& order) {
    ks_clear(o);
    switch (order) {
        case HtsSortOrder::UNSORTED:    kputs("unsorted", &o);   break;
        case HtsSortOrder::QUERYNAME:   kputs("queryname", &o);  break;
        case HtsSortOrder::COORDINATE:  kputs("coordinate", &o); break;
        default:                        kputs("unknown", &o);    break;
    }
    return o;
};
void operator>>(const char* s, HtsSortOrder& order) {
    if(s == NULL)                       order = HtsSortOrder::UNKNOWN;
    else if(!strcmp(s, "unsorted"))     order = HtsSortOrder::UNSORTED;
    else if(!strcmp(s, "queryname"))    order = HtsSortOrder::QUERYNAME;
    else if(!strcmp(s, "coordinate"))   order = HtsSortOrder::COORDINATE;
    else                                order = HtsSortOrder::UNKNOWN;
};

ostream& operator<<(ostream& o, const HtsGrouping& grouping) {
    switch (grouping) {
        case HtsGrouping::QUERY:        o << "query";       break;
        case HtsGrouping::REFERENCE:    o << "reference";   break;
        default:                        o << "none";
    }
    return o;
};
string& operator<<(string& o, const HtsGrouping& grouping) {
    switch (grouping) {
        case HtsGrouping::QUERY:        o.assign("query");      break;
        case HtsGrouping::REFERENCE:    o.assign("reference");  break;
        default:                        o.assign("none");       break;
    }
    return o;
};
kstring_t& operator<<(kstring_t& o, const HtsGrouping& grouping) {
    ks_clear(o);
    switch (grouping) {
        case HtsGrouping::QUERY:        kputs("query", &o);      break;
        case HtsGrouping::REFERENCE:    kputs("reference", &o);  break;
        default:                        kputs("unknown", &o);    break;
    }
    return o;
};
void operator>>(const char* s, HtsGrouping& grouping) {
    if(s == NULL)                       grouping = HtsGrouping::NONE;
    else if(!strcmp(s, "query"))        grouping = HtsGrouping::QUERY;
    else if(!strcmp(s, "reference"))    grouping = HtsGrouping::REFERENCE;
    else                                grouping = HtsGrouping::NONE;
};

ostream& operator<<(ostream& o, const Platform& platform) {
    switch (platform) {
        case Platform::CAPILLARY:   o << "CAPILLARY";   break;
        case Platform::LS454:       o << "LS454";       break;
        case Platform::ILLUMINA:    o << "ILLUMINA";    break;
        case Platform::SOLID:       o << "SOLID";       break;
        case Platform::HELICOS:     o << "HELICOS";     break;
        case Platform::IONTORRENT:  o << "IONTORRENT";  break;
        case Platform::ONT:         o << "ONT";         break;
        case Platform::PACBIO:      o << "PACBIO";      break;
        default:                    o << "UNKNOWN";     break;
    }
    return o;
};
string& operator<<(string& o, const Platform& platform) {
    switch (platform) {
        case Platform::CAPILLARY:   o.assign("CAPILLARY");  break;
        case Platform::LS454:       o.assign("LS454");      break;
        case Platform::ILLUMINA:    o.assign("ILLUMINA");   break;
        case Platform::SOLID:       o.assign("SOLID");      break;
        case Platform::HELICOS:     o.assign("HELICOS");    break;
        case Platform::IONTORRENT:  o.assign("IONTORRENT"); break;
        case Platform::ONT:         o.assign("ONT");        break;
        case Platform::PACBIO:      o.assign("PACBIO");     break;
        default:                    o.assign("UNKNOWN");    break;
    }
    return o;
};
kstring_t& operator<<(kstring_t& o, const Platform& platform) {
    ks_clear(o);
    switch (platform) {
        case Platform::CAPILLARY:   kputs("CAPILLARY", &o);  break;
        case Platform::LS454:       kputs("LS454", &o);      break;
        case Platform::ILLUMINA:    kputs("ILLUMINA", &o);   break;
        case Platform::SOLID:       kputs("SOLID", &o);      break;
        case Platform::HELICOS:     kputs("HELICOS", &o);    break;
        case Platform::IONTORRENT:  kputs("IONTORRENT", &o); break;
        case Platform::ONT:         kputs("ONT", &o);        break;
        case Platform::PACBIO:      kputs("PACBIO", &o);     break;
        default:                    kputs("UNKNOWN", &o);    break;
    }
    return o;
};
void operator>>(const char* s, Platform& platform) {
    if(s == NULL)                       platform = Platform::UNKNOWN;
    else if(!strcmp(s, "CAPILLARY"))    platform = Platform::CAPILLARY;
    else if(!strcmp(s, "LS454"))        platform = Platform::LS454;
    else if(!strcmp(s, "ILLUMINA"))     platform = Platform::ILLUMINA;
    else if(!strcmp(s, "SOLID"))        platform = Platform::SOLID;
    else if(!strcmp(s, "HELICOS"))      platform = Platform::HELICOS;
    else if(!strcmp(s, "IONTORRENT"))   platform = Platform::IONTORRENT;
    else if(!strcmp(s, "ONT"))          platform = Platform::ONT;
    else if(!strcmp(s, "PACBIO"))       platform = Platform::PACBIO;
    else                                platform = Platform::UNKNOWN;
};

ostream& operator<<(ostream& o, const FormatKind& kind) {
    switch (kind) {
        case FormatKind::FASTQ: o << "FASTQ";   break;
        case FormatKind::HTS:   o << "HTS";     break;
        default:                o << "UNKNOWN"; break;
    }
    return o;
};
string& operator<<(string& o, const FormatKind& kind) {
    switch (kind) {
        case FormatKind::FASTQ: o.assign("FASTQ");      break;
        case FormatKind::HTS:   o.assign("HTS");        break;
        default:                o.assign("UNKNOWN");    break;
    }
    return o;
};
void operator>>(const char* s, FormatKind& kind) {
    if(s == NULL)                   kind = FormatKind::UNKNOWN;
    else if(!strcmp(s, "FASTQ"))    kind = FormatKind::FASTQ;
    else if(!strcmp(s, "HTS"))      kind = FormatKind::HTS;
    else                            kind = FormatKind::UNKNOWN;
};

ostream& operator<<(ostream& o, const Decoder& decoder) {
    switch (decoder) {
        case Decoder::UNKNOWN:      o << "unknown";     break;
        case Decoder::MDD:          o << "mdd";         break;
        case Decoder::PAMLD:        o << "pamld";       break;
        case Decoder::BENCHMARK:    o << "benchmark";   break;
        default:                    o.setstate(ios_base::failbit);
    }
    return o;
};
string& operator<<(string& o, const Decoder& decoder) {
    switch (decoder) {
        case Decoder::UNKNOWN:      o.assign("unknown");    break;
        case Decoder::MDD:          o.assign("mdd");        break;
        case Decoder::PAMLD:        o.assign("pamld");      break;
        case Decoder::BENCHMARK:    o.assign("benchmark");  break;
        default:                                            break;
    }
    return o;
};
void operator>>(const char* s, Decoder& decoder) {
    if(s == NULL)                       decoder = Decoder::UNKNOWN;
    else if(!strcmp(s, "mdd"))          decoder = Decoder::MDD;
    else if(!strcmp(s, "pamld"))        decoder = Decoder::PAMLD;
    else if(!strcmp(s, "benchmark"))    decoder = Decoder::BENCHMARK;
    else                                decoder = Decoder::UNKNOWN;
};

ostream& operator<<(ostream& o, const LeftTokenOperator& operation) {
    switch (operation) {
        case LeftTokenOperator::NONE:               o << "none";                break;
        case LeftTokenOperator::REVERSE_COMPLEMENT: o << "reverse complement";  break;
    }
    return o;
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

ostream& operator<<(ostream& o, const htsExactFormat& hts_exact_format) {
    switch (hts_exact_format) {
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

ostream& operator<<(ostream& o, const htsCompression& hts_compression) {
    switch (hts_compression) {
        case htsCompression::gzip:      o << "gzip";            break;
        case htsCompression::bgzf:      o << "bgzf";            break;
        case htsCompression::custom:    o << "custom";          break;
        default:                        o << "no compression";  break;
    }
    return o;
};
