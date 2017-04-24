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

#ifndef PHENIQS_MODEL_H
#define PHENIQS_MODEL_H

#include <set>
#include <unordered_map>

#include <htslib/hts.h>
#include <htslib/hfile.h>
#include <htslib/kstring.h>
#include <rapidjson/document.h>

#include "constant.h"
#include "error.h"

using std::set;
using std::copy;
using std::hash;
using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::size_t;
using std::string;
using std::vector;
using std::ostream;
using std::ifstream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::make_pair;
using std::setprecision;
using std::unordered_map;

using rapidjson::Document;
using rapidjson::Value;
using rapidjson::SizeType;
using rapidjson::StringRef;

#define tag_to_code(t) uint16_t(*(t))<<8 | uint8_t(*((t) + 1))

/*  URL
*/
class URL {
    friend ostream& operator<<(ostream& o, const URL& url);
    friend string operator+(const string& lhs, const URL& rhs);
    friend string operator+(const URL& lhs, const string& rhs);
    friend bool operator<(const URL& lhs, const URL& rhs);

    public:
        URL();
        ~URL();
        URL(const URL& other);
        URL(const string& path, const IoDirection& direction);
        URL(const char* path, const size_t size, const IoDirection& direction);
        inline void clear() {
            ks_clear(_path);
            ks_clear(_name);
            ks_clear(_directory);
            ks_clear(_extension);
            ks_clear(_compression);
        };
        inline const char* const path() const {
            return _path.s;
        };
        inline const char* const name() const {
            return _name.s;
        };
        inline const char* const directory() const {
            return _directory.s;
        };
        inline const char* const extension() const {
            return _extension.s;
        };
        inline const char* const compression() const {
            return _compression.s;
        };
        inline const FormatType& type() const {
            return _type;
        };
        inline bool empty() const {
            return _path.l == 0;
        };
        inline bool is_file() const {
            return _name.l > 0;
        };
        inline bool is_directory() const {
            return _name.l == 0 && _directory.l > 0;
        };
        inline bool is_stdin() const {
            return !strcmp(_path.s, CANONICAL_STDIN_PATH);
        };
        inline bool is_stdout() const {
            return !strcmp(_path.s, CANONICAL_STDOUT_PATH);
        };
        inline bool is_stderr() const {
            return !strcmp(_path.s, CANONICAL_STDERR_PATH);
        };
        inline bool is_null() const {
            return !strcmp(_path.s, CANONICAL_NULL_DEVICE_PATH);
        };
        inline bool is_standard_stream() const {
            return is_stdin() || is_stdout() || is_stderr() || is_null();
        };
        inline bool is_absolute() const {
            return _directory.l > 0 && _directory.s[0] == PATH_SEPARATOR;
        };
        void parse(const char* path, const size_t size, const IoDirection& direction);
        void set_name(const char* name, const size_t size);
        void set_directory(const char* directory, const size_t size);
        void set_compression(const char* compression, const size_t size);
        void set_type(const char* type);
        void set_type(const FormatType type);
        void relocate(const URL& base);
        FormatKind kind() const;
        bool is_readable() const;
        bool is_writable() const;
        const char* const c_str() const;
        void describe(ostream& o) const;
        bool operator==(const URL& other) const;
        URL& operator=(const URL& other);
        operator string() const;
        void encode(Document& document, Value& value) const;

    private:
        kstring_t _path;
        kstring_t _name;
        kstring_t _directory;
        kstring_t _extension;
        kstring_t _compression;
        FormatType _type;

        inline void initialize() {
            ks_terminate(_path);
            ks_terminate(_name);
            ks_terminate(_directory);
            ks_terminate(_extension);
            ks_terminate(_compression);
        };
        void refresh();
        void expand(kstring_t* path);
        void decode_extension(const FormatType& type);
};
namespace std {
    template <> struct hash<URL> {
        size_t operator()(const URL& url) const {
            return hash<string>()(url.path());
        };
    };
};
/* Transform
*/
class Token {
    friend ostream& operator<<(ostream& o, const Token& token);
    void operator=(Token const &) = delete;

    public:
        const size_t index;
        const size_t input_segment_index;

        Token(
            const size_t& index,
            const size_t& input_segment_index,
            const int32_t& start,
            const int32_t& end,
            const bool& end_terminated);
        Token(const Token& other);
        string description() const;
        inline size_t decode_end(const size_t& length) const {
            size_t value;
            if(end_terminated) {
                int32_t v = (end < 0 ? length + end : end);
                value = v < 0 ? 0 : v;
                if(value > length) {
                    value = length;
                }
            } else {
                value = length;
            }
            return value;
        };
        inline size_t decode_start(const size_t& length) const {
            size_t value;
            int32_t v = start < 0 ? length + start : start;
            value = v < 0 ? 0 : v;
            if(value > length) {
                value = length;
            }
            return value;
        };
        inline bool empty() const {
            return (end_terminated && start >= end) && ((start >= 0 && end >= 0) || (start < 0 && end < 0));
        };
        inline bool constant() const {
            if(end_terminated) {
               return (start >= 0 && end >= 0) || (start < 0 && end < 0);
            } else {
                return start < 0;
            }
        };
        int32_t length() const {
            if(constant()) {
                if(end_terminated) {
                    return empty() ? 0 : end - start;
                } else {
                    return -start;
                }
            } else {
                return -1;
            }
        };
        operator string() const;

    private:
        const int32_t start;
        const int32_t end;
        const bool end_terminated;
};
class Transform {
    friend ostream& operator<<(ostream& o, const Transform& transform);
    void operator=(Transform const &) = delete;

    public:
        const size_t index;
        const size_t output_segment_index;
        const Token token;
        const LeftTokenOperator left;

        Transform(const size_t& index, const Token& token, const size_t& output_segment_index, const LeftTokenOperator& left);
        Transform(const Transform& other);
        string description() const;
        operator string() const;
};
/* DNA sequence
*/
class Sequence {
    friend ostream& operator<<(ostream& o, const Sequence& sequence);
    friend bool operator<(const Sequence& left, const Sequence& right);
    friend bool operator>(const Sequence& left, const Sequence& right);

    public:
        uint8_t* code;
        uint8_t* quality;
        size_t capacity;
        size_t length;

        Sequence();
        Sequence(const Sequence& other);
        Sequence& operator=(const Sequence& other);
        ~Sequence();
        void expected_error(float& error) const;
        void mask(const uint8_t& threshold);
        void fill(const char* code, const size_t& size);
        void fill(const uint8_t* code, const uint8_t* quality, const size_t& size);
        void append(const uint8_t* code, const uint8_t* quality, const size_t& size);
        void append(const Sequence& other, const size_t& start, const size_t& size);
        size_t append(const Sequence& other, const Transform& transform);
        inline void clear() {
            length = 0;
            code[length] = '\0';
            quality[length] = '\0';
        };
        inline size_t distance_from(const Sequence& other, const uint8_t threshold) const {
            size_t distance = 0;
            if(threshold > 0) {
                for (size_t i = 0; i < length; i++) {
                    if(other.quality[i] < threshold) {
                        // if quality is bellow threshold always count a miss
                        distance++;
                    } else if (code[i] != other.code[i]) {
                        distance++;
                    }
                }

            } else {
                for (size_t i = 0; i < length; i++) {
                    if (code[i] != other.code[i]) {
                        distance++;
                    }
                }
            }
            return distance;
        };
        inline size_t distance_from(const Sequence& other) const {
            size_t distance = 0;
            for (size_t i = 0; i < length; i++) {
                if (code[i] != other.code[i]) {
                    distance++;
                }
            }
            return distance;
        };
        inline string iupac_ambiguity() const {
            string result;
            result.reserve(length);
            for (size_t i = 0; i < length; i++) {
                result.push_back(BamToAmbiguousAscii[uint8_t(code[i])]);
            }
            return result;
        };
        inline void encode_iupac_ambiguity(kstring_t* buffer) const {
            if(length > 0) {
                ks_resize(buffer, buffer->l + length + 1);
                for (size_t i = 0; i < length; i++) {
                    buffer->s[buffer->l + i] = BamToAmbiguousAscii[uint8_t(code[i])];
                    // kputc(BamToAmbiguousAscii[uint8_t(code[i])], buffer);
                }
                buffer->l += length;
                ks_terminate(*buffer);
            }
        };
        inline void encode_iupac_ambiguity(string& buffer) const {
            if(length > 0) {
                buffer.reserve(buffer.size() + length + 1);
                for (size_t i = 0; i < length; i++) {
                    buffer.push_back(BamToAmbiguousAscii[uint8_t(code[i])]);
                }
            }
        };
        inline void encode_iupac_ambiguity(Value& value) const {
            char* buffer = (char*)malloc(length + 1);
            for (size_t i = 0; i < length; i++) {
                buffer[i] = BamToAmbiguousAscii[uint8_t(code[i])];
            }
            buffer[length] = '\0';
            value.SetString(StringRef(buffer, length));
        };
        inline void encode_phred_quality(kstring_t* buffer, const uint8_t phred_offset) const {
            if(length > 0) {
                ks_resize(buffer, buffer->l + length + 1);
                for (size_t i = 0; i < length; i++) {
                    buffer->s[buffer->l + i] = quality[i] + phred_offset;
                    // kputc(quality[i] + phred_offset, buffer);
                }
                buffer->l += length;
                ks_terminate(*buffer);
            }
        };
        inline void encode_phred_quality(string& buffer, const uint8_t phred_offset) const {
            if(length > 0) {
                buffer.reserve(buffer.size() + length + 1);
                for (size_t i = 0; i < length; i++) {
                    buffer.push_back(quality[i] + phred_offset);
                }
            }
        };
    private:
        inline void terminate() {
            code[length] = '\0';
            quality[length] = '\0';
        };
};
class Barcode {
    friend ostream& operator<<(ostream& o, const Barcode& barcode);

    public:
        Barcode();
        Barcode(const size_t& width);
        Barcode(const Barcode& other);
        Barcode& operator=(const Barcode& other);
        operator string() const;
        string iupac_ambiguity() const;
        string iupac_ambiguity(const size_t position) const;
        void append(const size_t& position, const Sequence& sequence, const Transform& transform);
        void fill(const size_t& position, const char* code, const size_t& size);
        void set_tolerance(const vector<uint8_t>& tolerance);
        void set_threshold(const uint8_t& threshold);
        inline void clear() {
            for(auto& fragment : fragments) {
                fragment.clear();
            }
            length = 0;
        };
        inline bool empty() const {
            return length == 0;
        };
        inline size_t size(const size_t& position) {
            return position < fragments.size() ? fragments[position].length : 0;
        };
        inline void resize(const size_t& width) {
            fragments.resize(width);
            tolerance.resize(width);
        };
        inline size_t total_fragments() const {
            return fragments.size();
        };
        inline bool corrected_match(const Barcode& other, size_t& distance) const {
            bool result = true;
            distance = 0;
            for (size_t i = 0; i < fragments.size(); i++) {
                size_t error = fragments[i].distance_from(other.fragments[i], threshold);
                if (error > tolerance[i]) {
                    result = false;
                }
                distance += error;
            }
            return result;
        };
        inline void encode_iupac_ambiguity(kstring_t* buffer) const {
            for(const auto& sequence : fragments) {
                sequence.encode_iupac_ambiguity(buffer);
            }
        };
        inline void encode_phred_quality(kstring_t* buffer, const uint8_t phred_offset) const {
            for(const auto& sequence : fragments) {
                sequence.encode_phred_quality(buffer, phred_offset);
            }
        };
        inline void decoding_probability(const Barcode& observed, double& probability, size_t& distance) const {
            double p = 1;
            size_t d = 0;
            for(size_t i = 0; i < fragments.size(); i++) {
                const Sequence& reference = fragments[i];
                const Sequence& o = observed.fragments[i];
                for(size_t j = 0; j < reference.length; j++) {
                    if (o.code[j] == reference.code[j]) {
                        p *= quality_to_inverse_probability(o.quality[j]);
                    } else {
                        d += 1;
                        if (o.code[j] != ANY_NUCLEOTIDE) {
                            p *= quality_to_third_probability(o.quality[j]);
                        } else {
                            p *= UNIFORM_BASE_PROBABILITY;
                        }
                    }
                }
            }
            distance = d;
            probability = p;
        };
        inline void accurate_decoding_probability(const Barcode& observed, double& probability, size_t& distance) const {
            double q = 0;
            size_t d = 0;
            for(size_t i = 0; i < fragments.size(); i++) {
                const Sequence& reference = fragments[i];
                const Sequence& o = observed.fragments[i];
                for(size_t j = 0; j < reference.length; j++) {
                    if (o.code[j] == reference.code[j]) {
                        q += quality_to_inverse_quality(o.quality[j]);
                    } else {
                        d += 1;
                        if (o.code[j] != ANY_NUCLEOTIDE) {
                            q += double(o.quality[j]);
                        } else {
                            q += UNIFORM_BASE_PHRED;
                        }
                    }
                }
            }
            distance = d;
            probability = pow(10.0, q * -0.1);
        };
        inline void compensated_decoding_probability(const Barcode& observed, double& probability, size_t& distance) const {
            // use the Kahan summation algorithm to minimize floating point drift
            // see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
            double sigma = 0;
            double compensation = 0;
            double y = 0;
            double t = 0;
            double q = 0;
            size_t d = 0;
            for(size_t i = 0; i < fragments.size(); i++) {
                const Sequence& reference = fragments[i];
                const Sequence& o = observed.fragments[i];
                for(size_t j = 0; j < reference.length; j++) {
                    if (o.code[j] == reference.code[j]) {
                        q = quality_to_inverse_quality(o.quality[j]);
                    } else {
                        d += 1;
                        if (o.code[j] != ANY_NUCLEOTIDE) {
                            q = double(o.quality[j]);
                        } else {
                            q = UNIFORM_BASE_PHRED;
                        }
                    }
                    y = q - compensation;
                    t = sigma + y;
                    compensation = (t - sigma) - y;
                    sigma = t;
                }
            }
            distance = d;
            probability = pow(10.0, sigma * -0.1);
        };
        void encode_configuration(Document& document, Value& node, const string& key) const;
        void encode_report(Document& document, Value& node, const string& key) const;
    private:
        size_t length;
        uint8_t threshold;
        vector< uint8_t > tolerance;
        vector< Sequence > fragments;
};
/*  @HD The header line

    VN  Format version
    SO  Sorting order of alignments
        unknown
        unsorted
        queryname
        coordinate
    GO  Grouping of alignments
        Indicating that similar alignment records are grouped together but the file is not necessarily sorted overall.
            none
            query       grouped by QNAME
            reference   grouped by RNAME/POS
*/
class HeadHDAtom {
    friend class HtsHeader;
    friend ostream& operator<<(ostream& o, const HeadHDAtom& hd);

    public:
        kstring_t VN;
        kstring_t SO;
        kstring_t GO;

        HeadHDAtom();
        ~HeadHDAtom();
        HeadHDAtom(const HeadHDAtom& other);
        HeadHDAtom& operator=(const HeadHDAtom& other);
        void set_alignment_sort_order(const HtsSortOrder& order);
        void set_alignment_grouping(const HtsGrouping& grouping);
        void set_version(const htsFormat* format);

    private:
        void encode(kstring_t* buffer) const;
        char* decode(char* position, const char* end);
};
/*  @PG Program

    ID  Program record identifier.
        Each @PG line must have a unique ID.
        The value of ID is used in the alignment PG tag and PP tags of other @PG lines.
        PG IDs may be modified when merging SAM files in order to handle collisions.
    PN  Program name
    CL  Command line
    PP  Previous @PG:ID. Must be valid reference to another @PG:ID
    DS  Description.
    VN  Program version
*/
class HeadPGAtom {
    friend class HtsHeader;
    friend ostream& operator<<(ostream& o, const HeadPGAtom& program);

    public:
        kstring_t ID;
        kstring_t PN;
        kstring_t CL;
        kstring_t PP;
        kstring_t DS;
        kstring_t VN;

        HeadPGAtom();
        ~HeadPGAtom();
        HeadPGAtom(const HeadPGAtom& other);
        HeadPGAtom& operator=(const HeadPGAtom& other);
        operator string() const;

    private:
        void encode(kstring_t* buffer) const;
        char* decode(char* position, const char* end);
};
/*  @RG Read Group

    ID  Read group identifier
        Each @RG line must have a unique ID.
        The value of ID is used in the RG tags of alignment records.
        Must be unique among all read groups in header section.
        Read group IDs may be modified when merging SAM files in order to handle collisions.
    LB  Library
    SM  Sample. Use pool name where a pool is being sequenced.
    PU  Platform unit Unique identifier. e.g. flowcell-barcode.lane for Illumina.
    CN  Name of sequencing center producing the read
    DS  Description
    DT  Date the run was produced ISO8601 date or datetime.
    PI  Predicted median insert size.
    PL  Platform or technology used to produce the reads.
            CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT, ONT, PACBIO
    PM  Platform model
    PG  Programs used for processing the read group
    FO  Flow order.
        The array of nucleotide bases that correspond to the nucleotides used for each flow of each read.
        Multi-base flows are encoded in IUPAC format, and non-nucleotide flows by various other characters.
        Format: /\*|[ACMGRSVTWYHKDBN]+/
    KS  The array of nucleotide bases that correspond to the key sequence of each read.

    From GATK: https://software.broadinstitute.org/gatk/guide/article?id=6472
    ID  Read group identifier This tag identifies which read group each read belongs to, so each read group's ID must be unique.
        It is referenced both in the read group definition line in the file header (starting with @RG) and in the RG:Z tag for each read record.
        Note that some Picard tools have the ability to modify IDs when merging SAM files in order to avoid collisions.
        In Illumina data, read group IDs are composed using the flowcell + lane name and number,
        making them a globally unique identifier across all sequencing data in the world.
        Use for BQSR: ID is the lowest denominator that differentiates factors contributing to technical batch effects:
        therefore, a read group is effectively treated as a separate run of the instrument in data processing steps
        such as base quality score recalibration, since they are assumed to share the same error model.

    PU  Platform Unit The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}.
        The {FLOWCELL_BARCODE} refers to the unique identifier for a particular flow cell.
        The {LANE} indicates the lane of the flow cell
        and the {SAMPLE_BARCODE} is a sample/library-specific identifier.
        Although the PU is not required by GATK but takes precedence over ID for base recalibration if it is present.
        In the example shown earlier, two read group fields, ID and PU, appropriately differentiate flow cell lane,
        marked by .2, a factor that contributes to batch effects.

    SM  Sample The name of the sample sequenced in this read group.
        GATK tools treat all read groups with the same SM value as containing sequencing data for the same sample,
        and this is also the name that will be used for the sample column in the VCF file.
        Therefore it's critical that the SM field be specified correctly.
        When sequencing pools of samples, use a pool name instead of an individual sample name.

    PL  Platform/technology used to produce the read This constitutes the only way to know what
        sequencing technology was used to generate the sequencing data. Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO.

    LB  DNA preparation library identifier MarkDuplicates uses the LB field to determine which read groups might contain molecular duplicates,
        in case the same DNA library was sequenced on multiple lanes.

    also informative: http://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
*/
class HeadRGAtom {
    friend class HtsHeader;
    friend ostream& operator<<(ostream& o, const HeadRGAtom& read_group);

    public:
        kstring_t ID;
        kstring_t PI;
        kstring_t LB;
        kstring_t SM;
        kstring_t PU;
        kstring_t CN;
        kstring_t DS;
        kstring_t DT;
        kstring_t PL;
        kstring_t PM;
        kstring_t PG;
        kstring_t FO;
        kstring_t KS;

        HeadRGAtom();
        ~HeadRGAtom();
        HeadRGAtom(const HeadRGAtom& other);
        HeadRGAtom& operator=(const HeadRGAtom& other);
        operator string() const;
        void set_platform(const Platform& value);
        void expand(const HeadRGAtom& other);
        void encode(Document& document, Value& node, const string& key) const;

    private:
        void encode(kstring_t* buffer) const;
        char* decode(char* position, const char* end);
};
/*  @CO free text comment
*/
class HeadCOAtom {
    friend class HtsHeader;
    friend ostream& operator<<(ostream& o, const HeadCOAtom& co);

    public:
        kstring_t CO;

        HeadCOAtom();
        ~HeadCOAtom();
        HeadCOAtom(const HeadCOAtom& other);
        HeadCOAtom& operator=(const HeadCOAtom& other);

    private:
        void encode(kstring_t* buffer) const;
        char* decode(char* position, const char* end);
};
/*  Feed specification
*/
class FeedSpecification {
    friend ostream& operator<<(ostream& o, const FeedSpecification& specification);

    public:
        const IoDirection direction;
        const size_t index;
        URL url;
        Platform platform;
        size_t capacity;
        size_t resolution;
        uint8_t phred_offset;
        unordered_map< string, const HeadPGAtom > program_by_id;
        unordered_map< string, const HeadRGAtom > read_group_by_id;
        hFILE* hfile;

        FeedSpecification (
            const IoDirection& direction,
            const size_t& index,
            const URL& url,
            const Platform& platform,
            const uint8_t& phred_offset);
        void set_capacity(const size_t& capacity);
        void set_resolution(const size_t& resolution);
        void register_rg(const HeadRGAtom& rg);
        void register_pg(const HeadPGAtom& pg);
        void describe(ostream& o) const;
        void probe();
};
/*  Channel specification
*/
class ChannelSpecification {
    friend ostream& operator<<(ostream& o, const ChannelSpecification& channel);

    public:
        size_t index;
        size_t TC;
        kstring_t FS;
        kstring_t CO;
        bool disable_quality_control;
        bool long_read;
        bool include_filtered;
        bool undetermined;
        double concentration;
        Barcode multiplex_barcode;
        vector< URL > output_urls;
        HeadRGAtom rg;

        ChannelSpecification(size_t index);
        ~ChannelSpecification();
        inline bool writable() const {
            return output_urls.size() > 0;
        };
        string alias() const;
        void describe(ostream& o) const;
        void encode(Document& document, Value& node) const;
};

/*  Barcode distance metric
*/
class Distance {
    public:
        Distance() :
            _min_word_length(0),
            _max_word_length(0),
            _min_distance(0),
            _max_distance(0),
            _shannon_bound(0),
            _padding(0),
            _spacing(1) {
        };
        void load() {
            _words.clear();
            _cumulative.clear();
            if (!_index.empty()) {

                _max_word_length = 0;
                _min_word_length = numeric_limits<size_t>::max();
                for(const auto& word : _index) {
                    _min_word_length = MIN(_min_word_length, word.size());
                    _max_word_length = MAX(_max_word_length, word.size());
                    _words.push_back(word);
                }

                _max_distance = 0;
                _min_distance = numeric_limits<size_t>::max();
                _matrix.resize(height());
                _cumulative.resize(height());
                for(size_t i = 0; i < height(); i++) {
                    const string& row = word(i);
                    _matrix[i].resize(height());
                    for(size_t j = 0; j < height(); j++) {
                        const string& column = word(j);
                        if(i == j) {
                            _matrix[i][j] = 0;

                        } else if(i < j) {
                            size_t distance = hamming_distance(row, column);
                            _min_distance = MIN(_min_distance, distance);
                            _max_distance = MAX(_max_distance, distance);
                            _matrix[i][j] = distance;
                            _cumulative[i] += distance;
                            _cumulative[j] += distance;

                        } else {
                            _matrix[i][j] = shannon_bound(row, column);
                        }
                    }
                }
                _shannon_bound = ((_min_distance - 1) / 2);

                for(size_t i = 0; i < _cumulative.size(); i++) {
                    _cumulative[i] /= (height() * 2);
                }
                // We want to know how many digits are in the biggest value to be able to align the matrix
                _padding = _spacing;
                size_t digit = _max_distance;
                do {
                    digit /= 10;
                    _padding++;
                } while (digit != 0);
            }
        };
        void add(const string& word) {
            _index.insert(word);
        };
        inline size_t width() const {
            return _max_word_length;
        };
        inline size_t height() const {
            return _words.size();
        };
        inline size_t value(size_t i, size_t j) const {
            return _matrix[i][j];
        };
        inline const string& word(size_t i) const {
            return _words[i];
        };
        inline size_t cumulative(size_t i) const {
            return _cumulative[i];
        };
        inline size_t minimum_distance() const {
            return _min_distance;
        };
        inline size_t shannon_bound() const {
            return _shannon_bound;
        };
        inline bool empty() const {
            return _words.empty();
        };
        void describe(ostream& o) const {
            o << std::left;
            if(!empty()) {
                for(size_t i = 0; i < height(); i++) {
                    o << '\t';
                    for(size_t j = 0; j < height(); j++) {
                        o << setw(_padding)<< value(i, j);
                    }
                    o << word(i) << ' ' << setw(_padding) << cumulative(i) << endl;
                }
                o << endl;
            }
        };

    private:
        size_t _min_word_length;
        size_t _max_word_length;
        size_t _min_distance;
        size_t _max_distance;
        size_t _shannon_bound;
        size_t _padding;
        size_t _spacing;
        set< string > _index;
        vector < string > _words;
        vector < size_t > _cumulative;
        vector< vector< size_t > > _matrix;

        inline size_t hamming_distance(const string& left, const string& right) const {
            size_t result = 0;
            for (size_t i = 0; i < left.length(); i++) {
                if (left[i] != right[i]) {
                    result++;
                }
            }
            return result;
        };
        inline size_t shannon_bound(const string& left, const string& right) const {
            size_t result = hamming_distance(left, right);
            result = ((result - 1) / 2);
            return result;
        };
};


#endif /* PHENIQS_MODEL_H */