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

#ifndef PHENIQS_PIPELINE_H
#define PHENIQS_PIPELINE_H

#include <string>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <htslib/thread_pool.h>

#include "error.h"
#include "json.h"
#include "url.h"
#include "nucleotide.h"
#include "phred.h"
#include "atom.h"
#include "accumulate.h"
#include "environment.h"
#include "feed.h"
#include "fastq.h"
#include "hts.h"

using std::set;
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
using std::setprecision;
using std::unordered_map;
using std::numeric_limits;

using std::mutex;
using std::condition_variable;
using std::unique_lock;
using std::thread;

enum class FormatKind : uint8_t {
    UNKNOWN,
    FASTQ,
    HTS,
};
void to_string(const FormatKind& value, string& result);
bool from_string(const char* value, FormatKind& result);
void to_kstring(const FormatKind& value, kstring_t& result);
bool from_string(const string& value, FormatKind& result);
ostream& operator<<(ostream& o, const FormatKind& value);
void encode_key_value(const string& key, const FormatKind& value, Value& container, Document& document);

static inline FormatKind format_kind_from_type(const FormatType& type) {
    switch(type) {
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
};

class Pivot;
class Channel;
class Pipeline;

/*  Pivot */

class Pivot {
Pivot(Pivot const &) = delete;
void operator=(Pivot const &) = delete;

public:
    thread pivot_thread;
    const int32_t index;
    const ProgramAction action;
    const Decoder decoder;
    bool determined;
    bool filtered;
    double multiplex_probability;
    double conditioned_multiplex_probability;
    int32_t multiplex_distance;
    Barcode multiplex_barcode;
    Barcode molecular_barcode;
    vector< Segment > input;
    vector< Segment > output;
    Channel* decoded_multiplex_channel;
    PivotAccumulator accumulator;
    unordered_map< string, ChannelAccumulator > channel_accumulator_by_barcode;
    unordered_map< string, ChannelAccumulator > channel_accumulator_by_read_group_id;
    Pivot(Pipeline& pipeline);
    void start();
    inline void clear();

private:
    Pipeline& pipeline;
    const bool disable_quality_control;
    const bool long_read;
    const Segment* leading_segment;
    void run();
    inline void validate();
    inline void decode_with_mdd();
    inline void decode_with_pamld();
    inline void transform();
    inline void decode_tagged_channel();
    inline void encode_auxiliary();
    inline void copy_auxiliary();
    inline void push();
    inline void increment();
    inline void encode_mmd_benchmark_auxiliary();
    inline void encode_pamld_benchmark_auxiliary();
};

/*  Channel */

class Channel {
public:
    mutex channel_mutex;
    const int32_t index;
    const string barcode_key;
    const string rg_key;
    const double concentration;
    const Barcode multiplex_barcode;
    const bool disable_quality_control;
    const bool long_read;
    const bool include_filtered;
    const bool undetermined;
    const bool writable;
    const HeadRGAtom rg;
    list< Feed* > output_feed_by_order;
    vector< Feed* > output_feed_by_segment;
    ChannelAccumulator channel_accumulator;
    Channel(const Pipeline& pipeline, const ChannelSpecification& specification);
    ~Channel();
    inline void increment(Pivot& pivot);
    void push(Pivot& pivot);
    void finalize(const PipelineAccumulator& pipeline_accumulator);
    void encode(Document& document, Value& value) const;
};

/*  Pipeline */

class Pipeline {
friend class Pivot;

public:
    Environment& environment;
    const Document& instruction;
    const ProgramAction program_action;
    const Decoder decoder;
    const Platform platform;
    const int32_t leading_segment_index;
    const bool disable_quality_control;
    const bool long_read;
    const bool include_filtered;
    const double adjusted_multiplex_noise_probability;
    const double random_multiplex_barcode_probability;
    const double multiplex_confidence;
    const int threads;
    const int32_t input_segment_cardinality;
    const int32_t output_segment_cardinality;
    const int32_t multiplex_segment_cardinality;
    const int32_t molecular_segment_cardinality;
    const int32_t concatenated_multiplex_barcode_length;
    const int32_t concatenated_molecular_barcode_length;
    InputSpecification input_specification;
    list< ChannelSpecification > channel_specification_array;
    Pipeline(Environment& environment);
    ~Pipeline();
    void execute();
    bool pull(Pivot& pivot);

private:
    bool end_of_input;
    htsThreadPool thread_pool;
    vector< Token > token_array;
    list< Transform > template_transform_array;
    list< Transform > multiplex_barcode_transform_array;
    list< Transform > molecular_barcode_transform_array;
    list< Feed* > input_feed_by_index;
    list< Feed* > output_feed_by_index;
    vector< Feed* > input_feed_by_segment;
    Channel* undetermined;
    list< Channel* > channel_array;
    unordered_map< string, Channel* > channel_by_barcode;
    unordered_map< string, Channel* > channel_by_read_group_id;
    vector< Pivot* > pivot_array;
    PivotAccumulator* input_accumulator;
    PipelineAccumulator* output_accumulator;
    void validate_url_accessibility();
    void load_input();
    void load_output();
    void load_pivot_array();
    void start();
    void stop();
    void finalize();
    void encode_report(ostream& o) const;
    void encode_quality_report(ostream& o) const;
    void encode_demultiplex_report(ostream& o) const;
    void encode_input_report(Document& document, Value& value) const;
    void encode_output_report(Document& document, Value& value) const;

    void probe(const URL& url);
    void calibrate(const URL& url);
    ChannelSpecification* load_channel_from_rg(const HeadRGAtom& rg);

};

#endif /* PHENIQS_PIPELINE_H */