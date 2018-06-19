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

#ifndef PHENIQS_ACCUMULATE_H
#define PHENIQS_ACCUMULATE_H

#include "include.h"

#include "url.h"
#include "barcode.h"
#include "read.h"

class NucleotideAccumulator;
class CycleAccumulator;
class AveragePhreadAccumulator;
class SegmentAccumulator;
class InputAccumulator;
class ChannelAccumulator;
class PipelineAccumulator;
class OutputAccumulator;

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

class SegmentAccumulator {
    void operator=(SegmentAccumulator const &) = delete;

    public:
        const int32_t index;
        const URL url;
        const uint8_t phred_offset;
        const Platform platform;
        const int32_t resolution;
        int32_t capacity;
        int32_t shortest;
        vector < uint64_t > nucleic_acid_count_by_code;
        AveragePhreadAccumulator average_phred;
        vector< CycleAccumulator > cycle_by_index;
        SegmentAccumulator(const Value& ontology);
        SegmentAccumulator(const SegmentAccumulator& other) :
            index(other.index),
            url(other.url),
            phred_offset(other.phred_offset),
            platform(other.platform),
            resolution(other.resolution),
            capacity(other.capacity),
            shortest(other.shortest),
            nucleic_acid_count_by_code(other.nucleic_acid_count_by_code) {
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

class ChannelAccumulator {
    public:
        const uint64_t index;
        const Algorithm algorithm;
        const bool disable_quality_control;
        const bool undetermined;
        const double concentration;
        const Barcode barcode;
        const HeadRGAtom rg;
        uint64_t count;
        double multiplex_distance;
        double multiplex_confidence;
        uint64_t pf_count;
        double pf_multiplex_distance;
        double pf_multiplex_confidence;
        double pf_fraction;
        double pooled_fraction;
        double pf_pooled_fraction;
        double pooled_multiplex_fraction;
        double pf_pooled_multiplex_fraction;
        uint64_t accumulated_multiplex_distance;
        double accumulated_multiplex_confidence;
        uint64_t accumulated_pf_multiplex_distance;
        double accumulated_pf_multiplex_confidence;
        vector< SegmentAccumulator > segment_by_index;

        ChannelAccumulator(const Value& ontology);
        inline void increment(const Read& read) {
            ++count;
            if(read.multiplex_distance()) {
                accumulated_multiplex_distance += static_cast< uint64_t >(read.multiplex_distance());
            }
            if(algorithm == Algorithm::PAMLD) {
                accumulated_multiplex_confidence += read.multiplex_error();
            }
            if(!read.qcfail()) {
                ++pf_count;
                if(read.multiplex_distance()) {
                    accumulated_pf_multiplex_distance += static_cast< uint64_t >(read.multiplex_distance());
                }
                if(algorithm == Algorithm::PAMLD) {
                    accumulated_pf_multiplex_confidence += read.multiplex_error();
                }
            }
            for(size_t i(0); i < segment_by_index.size(); ++i) {
                segment_by_index[i].increment(read[i]);
            }
        };
        void finalize(const OutputAccumulator& decoder_accumulator);
        ChannelAccumulator& operator+=(const ChannelAccumulator& rhs);
};
template<> vector< ChannelAccumulator > decode_value_by_key(const Value::Ch* key, const Value& container);
bool encode_value(const ChannelAccumulator& value, Value& container, Document& document);
bool encode_key_value(const string& key, const ChannelAccumulator& value, Value& container, Document& document);

class InputAccumulator {
    public:
        const bool disable_quality_control;
        uint64_t count;
        uint64_t pf_count;
        double pf_fraction;
        vector< SegmentAccumulator > segment_by_index;
        InputAccumulator(const Value& ontology);
        inline void increment(const Read& read) {
            ++count;
            if(!read.qcfail()) {
                ++pf_count;
            }
            for(size_t i(0); i < segment_by_index.size(); ++i) {
                segment_by_index[i].increment(read[i]);
            }
        };
        void finalize();
        InputAccumulator& operator+=(const InputAccumulator& rhs);
};
bool encode_key_value(const string& key, const InputAccumulator& value, Value& container, Document& document);

class OutputAccumulator {
    public:
        const Algorithm algorithm;
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
        double accumulated_multiplex_confidence;
        uint64_t accumulated_pf_multiplex_distance;
        double accumulated_pf_multiplex_confidence;
        vector< ChannelAccumulator > channel_by_index;
        ChannelAccumulator undetermined;

        OutputAccumulator(const Value& ontology);
        inline void increment(const size_t& index, const Read& read) {
            if(index > 0) {
                channel_by_index[index - 1].increment(read);
            } else {
                undetermined.increment(read);
            }
        };
        void finalize();
        OutputAccumulator& operator+=(const OutputAccumulator& rhs);
};
bool encode_key_value(const string& key, const OutputAccumulator& value, Value& container, Document& document);

#endif /* PHENIQS_ACCUMULATE_H */
