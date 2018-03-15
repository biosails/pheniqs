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
#include "constant.h"
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

class Pivot;
class Channel;
class Pipeline;

/*  Pivot */

class Pivot {
Pivot(Pivot const &) = delete;
void operator=(Pivot const &) = delete;

public:
    thread pivot_thread;
    const size_t index;
    const ProgramAction action;
    const Decoder decoder;
    bool determined;
    bool filtered;
    double multiplex_probability;
    double conditioned_multiplex_probability;
    size_t multiplex_distance;
    Barcode multiplex_barcode;
    Barcode molecular_barcode;
    vector< Segment > input;
    vector< Segment > output;
    Channel* decoded_multiplex_channel;
    PivotAccumulator accumulator;
    unordered_map< string, ChannelAccumulator > channel_by_barcode;
    unordered_map< string, ChannelAccumulator > channel_by_read_group_id;
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
    const size_t index;
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
    vector< Feed* > output_feeds;
    vector< Feed* > unique_output_feeds;
    ChannelAccumulator channel_accumulator;
    Channel(const Pipeline& pipeline, const ChannelSpecification& specification);
    ~Channel();
    inline void increment(Pivot& pivot);
    void push(Pivot& pivot);
    void finalize(const uint64_t& pool_count, const uint64_t& pool_pf_count);
    void encode(Document& document, Value& value, const bool disable_quality_control) const;
};

/*  Pipeline */

class Pipeline {
friend class Pivot;

public:
    Environment& environment;
    Pipeline(Environment& environment);
    ~Pipeline();
    void execute();
    bool pull(Pivot& pivot);
    void finalize();
    void encode_report(ostream& o) const;
    Feed* resolve_feed(const URL& url, const IoDirection& direction) const;

private:
    const double conditioned_probability_threshold;
    const bool disable_quality_control;
    const bool long_read;
    bool end_of_input;
    htsThreadPool thread_pool;

    vector< Feed* > input_feeds;
    vector< Feed* > unique_input_feeds;
    vector< Feed* > unique_output_feeds;
    unordered_map< URL, Feed* > input_feed_by_url;
    unordered_map< URL, Feed* > output_feed_by_url;

    Channel* undetermined;
    vector< Channel* > channels;
    unordered_map< string, Channel* > channel_by_barcode;
    unordered_map< string, Channel* > channel_by_read_group_id;

    vector< Pivot* > pivots;
    PivotAccumulator* input_accumulator;
    PipelineAccumulator* output_accumulator;

    void start();
    void stop();
    void probe(const URL& url);
    Feed* load_feed(FeedSpecification* specification);
    void load_input_feeds();
    void load_output_feeds();
    inline void initialize_channels();
    inline void initialize_pivots();
    void encode_quality_report(ostream& o) const;
    void encode_demultiplex_report(ostream& o) const;
    void encode_input_report(Document& document, Value& value) const;
    void encode_output_report(Document& document, Value& value) const;
};

#endif /* PHENIQS_PIPELINE_H */