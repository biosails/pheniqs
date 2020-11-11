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

#ifndef PHENIQS_CHANNEL_H
#define PHENIQS_CHANNEL_H

#include "include.h"
#include "feed.h"

class AveragePhreadAccumulator {
    public:
        uint64_t count;
        double min_value;
        double max_value;
        double sum_value;
        double mean_value;
        vector< uint64_t > distribution;
        AveragePhreadAccumulator();
        AveragePhreadAccumulator(const AveragePhreadAccumulator& other) :
            count(other.count),
            min_value(other.min_value),
            max_value(other.max_value),
            sum_value(other.sum_value),
            mean_value(other.mean_value),
            distribution(other.distribution) {
        };
        inline void increment(const Segment& segment) {
            ++count;
            double value(0);
            for(int32_t i(0); i < segment.length; ++i) {
                value += segment.quality[i];
            }
            value /= double(segment.length);
            sum_value += value;
            min_value = min(min_value, value);
            max_value = max(max_value, value);
            ++(distribution[static_cast< size_t >(value)]);
        };
        void finalize();
        AveragePhreadAccumulator& operator=(const AveragePhreadAccumulator& rhs);
        AveragePhreadAccumulator& operator+=(const AveragePhreadAccumulator& rhs);
};

class NucleotideAccumulator {
    public:
        uint64_t count;
        uint8_t min_quality;
        uint8_t max_quality;
        uint64_t sum_quality;
        double mean_quality;
        uint8_t Q1;
        uint8_t Q3;
        uint8_t IQR;
        uint8_t LW;
        uint8_t RW;
        uint8_t median_quality;
        vector< uint64_t > distribution;
        NucleotideAccumulator();
        NucleotideAccumulator(const NucleotideAccumulator& other) :
            count(other.count),
            min_quality(other.min_quality),
            max_quality(other.max_quality),
            sum_quality(other.sum_quality),
            mean_quality(other.mean_quality),
            Q1(other.Q1),
            Q3(other.Q3),
            IQR(other.IQR),
            LW(other.LW),
            RW(other.RW),
            median_quality(other.median_quality),
            distribution(other.distribution) {
        };
        inline void increment(const uint8_t phred) {
            ++(distribution[phred]);
        };
        inline uint64_t quantile(const double portion) {
            uint64_t position(portion * count);
            uint8_t phred(0);
            while (position > 0) {
                if(distribution[phred] >= position) {
                    break;
                }
                position -= distribution[phred];
                ++phred;
                while (distribution[phred] == 0) {
                    ++phred;
                }
            }
            return phred;
        };
        void finalize();
        NucleotideAccumulator& operator=(const NucleotideAccumulator& rhs);
        NucleotideAccumulator& operator+=(const NucleotideAccumulator& rhs);
};

class CycleAccumulator {
    public:
        vector< NucleotideAccumulator > nucleotide_by_code;
        CycleAccumulator();
        CycleAccumulator(const CycleAccumulator& other) :
            nucleotide_by_code(other.nucleotide_by_code) {
        };
        inline void increment(const uint8_t nucleotide, const uint8_t phred) {
            nucleotide_by_code[nucleotide].increment(phred);
        };
        void finalize();
        CycleAccumulator& operator=(const CycleAccumulator& rhs);
        CycleAccumulator& operator+=(const CycleAccumulator& rhs);
};

class SegmentAccumulator {
    public:
        void operator=(SegmentAccumulator const &) = delete;
        int32_t capacity;
        int32_t shortest;
        vector < uint64_t > nucleic_acid_count_by_code;
        AveragePhreadAccumulator average_phred;
        vector< CycleAccumulator > cycle_by_index;
        SegmentAccumulator();
        SegmentAccumulator(const SegmentAccumulator& other) :
            // index(other.index),
            capacity(other.capacity),
            shortest(other.shortest),
            nucleic_acid_count_by_code(other.nucleic_acid_count_by_code),
            average_phred(other.average_phred),
            cycle_by_index(other.cycle_by_index) {
        };
        inline void increment(const Segment& segment) {
            if(segment.length > capacity) {
                cycle_by_index.resize(segment.length);
                capacity = segment.length;
            }
            if(segment.length < shortest) {
                shortest = segment.length;
            }
            for(int32_t i(0); i < segment.length; ++i) {
                ++(nucleic_acid_count_by_code[NO_NUCLEOTIDE]);
                ++(nucleic_acid_count_by_code[segment.code[i]]);
                cycle_by_index[i].increment(segment.code[i], segment.quality[i]);
            }
            average_phred.increment(segment);
        };
        void finalize();
        SegmentAccumulator& operator+=(const SegmentAccumulator& rhs);
};
bool encode_value(const SegmentAccumulator& value, Value& container, Document& document);

class ReadAccumulator {
    public:
        void operator=(ReadAccumulator const &) = delete;
        vector< SegmentAccumulator > segment_accumulator_by_index;
        ReadAccumulator(const int32_t& cardinality);
        ReadAccumulator(const ReadAccumulator& other) :
            segment_accumulator_by_index(other.segment_accumulator_by_index) {
        };
        inline void increment(const Read& read) {
            for(size_t i(0); i < segment_accumulator_by_index.size(); ++i) {
                segment_accumulator_by_index[i].increment(read[i]);
            }
        };
        void finalize();
        ReadAccumulator& operator+=(const ReadAccumulator& rhs);
};
bool encode_value(const ReadAccumulator& value, Value& container, Document& document);

class Channel {
    public:
        void operator=(Channel const &) = delete;
        const int32_t index;
        // const HeadRGAtom rg;
        const bool filter_outgoing_qc_fail;
        const bool enable_quality_control;
        const list< URL > output_feed_url_by_segment;
        vector< Feed* > output_feed_lock_order;
        vector< Feed* > output_feed_by_segment;
        ReadAccumulator read_accumulator;

        Channel(const Value& ontology);
        Channel(const Channel& other);
        inline void push(const Read& read) {
            if(output_feed_lock_order.size() > 0) {
                if(!filter_outgoing_qc_fail || !read.qcfail()) {
                    /* acquire a push lock for all feeds in a fixed order */
                    vector< unique_lock< mutex > > feed_locks;
                    feed_locks.reserve(output_feed_lock_order.size());
                    for(const auto feed : output_feed_lock_order) {
                        feed_locks.push_back(feed->acquire_push_lock());
                    }

                    /* push the segments to the output feeds */
                    for(size_t i(0); i < output_feed_by_segment.size(); ++i) {
                        output_feed_by_segment[i]->push(read[i]);
                    }

                    /* release the locks on the feeds in reverse order */
                    for(auto feed_lock(feed_locks.rbegin()); feed_lock != feed_locks.rend(); ++feed_lock) {
                        feed_lock->unlock();
                    }
                }
            }
            if(enable_quality_control) {
                read_accumulator.increment(read);
            }
        };
        void populate(unordered_map< URL, Feed* >& output_feed_by_url);
        void finalize();
        void encode(Value& container, Document& document) const;
        Channel& operator+=(const Channel& rhs);
};
template<> vector< Channel > decode_value< vector< Channel > >(const Value& container);

class Multiplexer {
    public:
        const bool filter_outgoing_qc_fail;
        const bool enable_quality_control;
        vector< Channel > channel_by_index;

        Multiplexer(const Value& ontology);
        Multiplexer(const Multiplexer& other);
        inline void push(const Read& read) {
            channel_by_index[read.channel_index].push(read);
        };
        void populate(unordered_map< URL, Feed* >& output_feed_by_url) {
            for(auto& channel : channel_by_index) {
                channel.populate(output_feed_by_url);
            }
        };
        void collect(const Multiplexer& other);
        void finalize();
        void encode(Value& container, Document& document) const;
};

#endif /* PHENIQS_CHANNEL_H */
