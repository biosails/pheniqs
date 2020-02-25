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
#include "json.h"
#include "phred.h"

class AccumulatingTag;
class AccumulatingClassifier;
class NucleotideAccumulator;
class CycleAccumulator;

class AccumulatingTag {
    public:
        uint64_t count;
        uint64_t pf_count;
        uint64_t accumulated_distance;
        double accumulated_confidence;
        uint64_t low_conditional_confidence_count;
        uint64_t low_confidence_count;
        uint64_t accumulated_pf_distance;
        double accumulated_pf_confidence;

        double pf_fraction;                     /*  pf_count / count */
        double average_distance;                /*  accumulated_distance / count */
        double average_confidence;              /*  accumulated_confidence / count */
        double average_pf_distance;             /*  accumulated_pf_distance / pf_count */
        double average_pf_confidence;           /*  accumulated_pf_confidence / pf_count */
        double pooled_fraction;                 /*  count / decoder.count */
        double pf_pooled_fraction;              /*  pf_count / decoder.pf_count */
        double pooled_classified_fraction;      /*  count / decoder.classified_count */
        double pf_pooled_classified_fraction;   /*  pf_count / decoder.pf_classified_count */

        AccumulatingTag();
        AccumulatingTag(const AccumulatingTag& other);
        virtual ~AccumulatingTag() {};
        virtual void finalize(const AccumulatingClassifier& parent);
        virtual void encode(Value& container, Document& document) const;
        AccumulatingTag& operator+=(const AccumulatingTag& rhs);
};

class AccumulatingClassifier {
    public:
        const int32_t index;
        uint64_t count;
        uint64_t pf_count;
        uint64_t classified_count;
        uint64_t accumulated_classified_distance;
        double accumulated_classified_confidence;
        uint64_t low_conditional_confidence_count;
        uint64_t low_confidence_count;
        uint64_t pf_classified_count;
        uint64_t accumulated_pf_classified_distance;
        double accumulated_pf_classified_confidence;

        double pf_fraction;                         /*  pf_count / count */
        double classified_fraction;                 /*  classified_count / count */
        double average_classified_distance;         /*  accumulated_classified_distance / classified_count */
        double average_classified_confidence;       /*  accumulated_classified_confidence / classified_count */
        double pf_classified_fraction;              /*  pf_classified_count / pf_count */
        double classified_pf_fraction;              /*  pf_classified_count / classified_count */
        double average_pf_classified_distance;      /*  accumulated_pf_classified_distance / pf_classified_count */
        double average_pf_classified_confidence;    /*  accumulated_pf_classified_confidence / pf_classified_count */

        AccumulatingClassifier(const int32_t index);
        AccumulatingClassifier(const AccumulatingClassifier& other);
        virtual ~AccumulatingClassifier() {};
        virtual void finalize();
        virtual void encode(Value& container, Document& document) const;
        AccumulatingClassifier& operator+=(const AccumulatingClassifier& rhs);
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

#endif /* PHENIQS_ACCUMULATE_H */
