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

#ifndef PHENIQS_ACCUMULATE_H
#define PHENIQS_ACCUMULATE_H

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
#include "nucleotide.h"
#include "phred.h"
#include "sequence.h"
#include "segment.h"
#include "model.h"

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
        size_t quantile(const double portion);
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
        uint64_t pf_count;
        double pf_fraction;
        vector< FeedAccumulator > feed_accumulators;
        PivotAccumulator(const size_t total_input_segments);
        inline void increment(const bool& filtered, const vector< Segment >& input) {
            count++;
            if (!filtered) {
                pf_count++;
            }
            for(size_t i = 0; i < feed_accumulators.size(); i++) {
                feed_accumulators[i].increment(input[i].sequence);
            }
        };
        void encode(Document& document, Value& value, const bool disable_quality_control) const;
        void finalize();
        PivotAccumulator& operator+=(const PivotAccumulator& rhs);
};
class ChannelAccumulator {
    public:
        const bool undetermined;
        uint64_t count;
        double multiplex_distance;
        double multiplex_confidence;
        uint64_t pf_count;
        double pf_multiplex_distance;
        double pf_multiplex_confidence;
        double pf_fraction;
        double pooled_fraction;
        double pf_pooled_fraction;
        uint64_t accumulated_multiplex_distance;
        uint64_t accumulated_multiplex_confidence;
        uint64_t accumulated_pf_multiplex_distance;
        uint64_t accumulated_pf_multiplex_confidence;
        vector< FeedAccumulator > feed_accumulators;
        ChannelAccumulator(const ChannelSpecification& specification);
        inline void increment (
            const Decoder& decoder,
            const bool& filtered,
            const double& multiplex_probability,
            const size_t& multiplex_distance,
            const vector< Segment >& output) {

            count++;
            if(multiplex_distance) {
                accumulated_multiplex_distance += multiplex_distance;
            }
            if(decoder == Decoder::PAMLD) {
                accumulated_multiplex_confidence += multiplex_probability; 
            }
            if (!filtered) {
                pf_count++;
                if(multiplex_distance) {
                    accumulated_pf_multiplex_distance += multiplex_distance;
                }
                if(decoder == Decoder::PAMLD) {
                    accumulated_pf_multiplex_confidence += multiplex_probability;
                }
            }

            for(size_t i = 0; i < feed_accumulators.size(); i++) {
                feed_accumulators[i].increment(output[i].sequence);
            }
        };
        void encode(Document& document, Value& value, const bool disable_quality_control) const;
        void finalize(const uint64_t& pool_count, const uint64_t& pool_pf_count);
        ChannelAccumulator& operator+=(const ChannelAccumulator& rhs);
};
class PipelineAccumulator {
    public:
        uint64_t count;
        uint64_t multiplex_count;
        double multiplex_fraction;
        double multiplex_distance;
        double multiplex_confidence;
        uint64_t pf_count;
        double pf_fraction;
        uint64_t pf_multiplex_count;
        double pf_multiplex_fraction;
        double pf_multiplex_distance;
        double pf_multiplex_confidence;
        double multiplex_pf_fraction;
        uint64_t accumulated_multiplex_distance;
        uint64_t accumulated_multiplex_confidence;
        uint64_t accumulated_pf_multiplex_distance;
        uint64_t accumulated_pf_multiplex_confidence;

        PipelineAccumulator();
        void collect(const ChannelAccumulator& channel_accumulator);
        void encode(Document& document, Value& value) const;
        void finalize();
};

#endif /* PHENIQS_ACCUMULATE_H */