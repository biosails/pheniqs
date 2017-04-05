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

#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/error/en.h>

#include "feed.h"
#include "environment.h"

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

using rapidjson::Document;
using rapidjson::Value;
using rapidjson::SizeType;
using rapidjson::StringBuffer;
using rapidjson::PrettyWriter;


class Pivot;
class Channel;
class Pipeline;

/*  Quality */

class NucleicAcidAccumulator {
    NucleicAcidAccumulator(NucleicAcidAccumulator const &) = delete;
    void operator=(NucleicAcidAccumulator const &) = delete;

    public:
        uint64_t count;
        uint64_t min;
        uint64_t max;
        uint64_t sum;
        double mean;
        uint64_t Q1;
        uint64_t Q3;
        uint64_t IQR;
        uint64_t LW;
        uint64_t RW;
        uint64_t median;
        uint64_t distribution[EFFECTIVE_PHRED_RANGE];
        NucleicAcidAccumulator();
        inline void increment(const uint8_t phred) {
            distribution[phred]++;
        };
        inline size_t quantile(const double portion);
        void finalize();
        NucleicAcidAccumulator& operator+=(const NucleicAcidAccumulator& rhs);
};
class CycleAccumulator {
    CycleAccumulator(CycleAccumulator const &) = delete;
    void operator=(CycleAccumulator const &) = delete;

    public:
        vector<NucleicAcidAccumulator> iupac_nucleic_acid;
        CycleAccumulator();
        inline void increment(const uint8_t nucleotide, const uint8_t phred) {
            iupac_nucleic_acid[nucleotide].increment(phred);
        };
        void finalize();
        CycleAccumulator& operator+=(const CycleAccumulator& rhs);
};
class SegmentAccumulator {
    SegmentAccumulator(SegmentAccumulator const &) = delete;
    void operator=(SegmentAccumulator const &) = delete;

    public:
        uint64_t count;
        double min;
        double max;
        double sum;
        double mean;
        uint64_t distribution[EFFECTIVE_PHRED_RANGE];
        SegmentAccumulator();
        inline void increment(const Sequence& sequence) {
            count++;
            double value = 0.0;
            for (size_t i = 0; i < sequence.length; i++) {
                value += sequence.quality[i];
            }
            value /= double(sequence.length);
            sum += value;
            min = MIN(min, value);
            max = MAX(max, value);
            distribution[(uint8_t)value]++;
        };
        void finalize();
        SegmentAccumulator& operator+=(const SegmentAccumulator& rhs);
};
class FeedAccumulator {
    FeedAccumulator(FeedAccumulator const &) = delete;
    void operator=(FeedAccumulator const &) = delete;

    public:
        uint64_t length;
        uint64_t shortest;
        uint64_t iupac_nucleic_acid_count[IUPAC_CODE_SIZE];
        SegmentAccumulator average_phred;
        vector< CycleAccumulator* > cycles;
        FeedAccumulator();
        ~FeedAccumulator();
        inline void increment(const Sequence& sequence) {
            if(sequence.length > length) {
                for(size_t i = length; i < sequence.length; i++) {
                    cycles.push_back(new CycleAccumulator());
                }
                length = sequence.length;
            }
            if(sequence.length < shortest) {
                shortest = sequence.length;
            }
            for(size_t i = 0; i < sequence.length; i++) {
                iupac_nucleic_acid_count[NO_NUCLEOTIDE]++;
                iupac_nucleic_acid_count[sequence.code[i]]++;
                cycles[i]->increment(sequence.code[i], sequence.quality[i]);
            }
            average_phred.increment(sequence);
        };
        void encode(Document& document, Value& value) const;
        void finalize();
        FeedAccumulator& operator+=(const FeedAccumulator& rhs);
};
class PivotAccumulator {
    public:
        uint64_t count;
        uint64_t filtered_count;
        uint64_t determined_count;
        uint64_t determined_filtered_count;
        vector< FeedAccumulator > feed_accumulators;
        PivotAccumulator(const size_t total_input_segments);
        inline void increment(const Pivot& pivot);
        void encode(Document& document, Value& value, const bool disable_quality_control) const;
        void finalize();
        PivotAccumulator& operator+=(const PivotAccumulator& rhs);
};
class ChannelAccumulator {
    public:
        uint64_t count;
        uint64_t filtered_count;
        uint64_t perfect_count;
        uint64_t yield;
        uint64_t perfect_yield;
        uint64_t filtered_yield;
        uint64_t distance_sum;
        double probability_sum;
        const bool include_filtered;
        vector< FeedAccumulator > feed_accumulators;
        ChannelAccumulator(const ChannelSpecification& specification);
        inline void increment(const Pivot& pivot);
        void encode(Document& document, Value& value, const bool disable_quality_control) const;
        void finalize();
        ChannelAccumulator& operator+=(const ChannelAccumulator& rhs);
};

/*  Pivot */

class Pivot {
    Pivot(Pivot const &) = delete;
    void operator=(Pivot const &) = delete;

    public:
        thread pivot_thread;
        const size_t index;
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
        Pivot(Pipeline& pipeline);
        void start();
        inline void clear();

    private:
        Pipeline& pipeline;
        const bool disable_quality_control;
        const bool long_read;
        const Segment* leading_segment;
        void run();
        void rotate();
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
        const string key;
        const double concentration;
        const Barcode multiplex_barcode;
        const bool disable_quality_control;
        const bool long_read;
        const bool include_filtered;
        const bool undetermined;
        const HeadRGAtom rg;
        vector< Feed* > output_feeds;
        vector< Feed* > unique_output_feeds;
        ChannelAccumulator channel_accumulator;
        Channel(const Pipeline& pipeline, const ChannelSpecification& specification);
        ~Channel();
        inline void increment(Pivot& pivot);
        void push(Pivot& pivot);
        void finalize();
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

        unordered_map< URL, Feed* > input_feed_by_url;
        unordered_map< URL, Feed* > output_feed_by_url;
        vector< Feed* > unique_input_feeds;
        vector< Feed* > unique_output_feeds;

        // input feeds
        vector< Feed* > input_feeds;

        // output channels
        Channel* undetermined;
        vector< Channel* > channels;
        unordered_map< string, Channel* > channel_by_barcode;
        unordered_map< string, Channel* > channel_by_read_group_id;

        // transforming processors
        vector< Pivot* > pivots;

        // accumulator to aggregate all pivot accumulators
        PivotAccumulator input_accumulator;

        void start();
        void stop();
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