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

#include "accumulator.h"

/*  AccumulatingTag */

AccumulatingTag::AccumulatingTag() :
    count(0),
    pf_count(0),
    accumulated_distance(0),
    accumulated_confidence(0),
    low_conditional_confidence_count(0),
    low_confidence_count(0),
    accumulated_pf_distance(0),
    accumulated_pf_confidence(0),

    pf_fraction(0),
    average_distance(0),
    average_confidence(0),
    average_pf_distance(0),
    average_pf_confidence(0),
    pooled_fraction(0),
    pf_pooled_fraction(0),
    pooled_classified_fraction(0),
    pf_pooled_classified_fraction(0) {
};
AccumulatingTag::AccumulatingTag(const AccumulatingTag& other) :
    count(other.count),
    pf_count(other.pf_count),
    accumulated_distance(other.accumulated_distance),
    accumulated_confidence(other.accumulated_confidence),
    low_conditional_confidence_count(other.low_conditional_confidence_count),
    low_confidence_count(other.low_confidence_count),
    accumulated_pf_distance(other.accumulated_pf_distance),
    accumulated_pf_confidence(other.accumulated_pf_confidence) {
};
void AccumulatingTag::finalize(const AccumulatingClassifier& parent) {
    if(count > 0) {
        average_distance = accumulated_distance / double(count);
        average_confidence = accumulated_confidence / double(count);
        if(parent.count > 0) {
            pooled_fraction = double(count) / double(parent.count);
        }
        if(parent.classified_count > 0) {
            pooled_classified_fraction = double(count) / double(parent.classified_count);
        }
    }
    if(pf_count > 0) {
        pf_fraction = double(pf_count) / double(count);
        average_pf_distance = accumulated_pf_distance / double(pf_count);
        average_pf_confidence = accumulated_pf_confidence / double(pf_count);

        if(parent.pf_count > 0) {
            pf_pooled_fraction = double(pf_count) / double(parent.pf_count);
        }
        if(parent.pf_classified_count > 0) {
            pf_pooled_classified_fraction = double(pf_count) / double(parent.pf_classified_count);
        }
    }
};
void AccumulatingTag::encode(Value& container, Document& document) const {
    encode_key_value("count", count, container, document);
    if(average_distance > 0) {
        encode_key_value("average distance", average_distance, container, document);
    }
    if(average_confidence > 0) {
        encode_key_value("average confidence", average_confidence, container, document);
    }
    if(low_conditional_confidence_count > 0) {
        encode_key_value("low conditional confidence count", low_conditional_confidence_count, container, document);
    }
    if(low_confidence_count > 0) {
        encode_key_value("low confidence count", low_confidence_count, container, document);
    }
    encode_key_value("pooled fraction", pooled_fraction, container, document);
    if(pooled_classified_fraction > 0) {
        encode_key_value("pooled classified fraction", pooled_classified_fraction, container, document);
    }
    encode_key_value("pf count", pf_count, container, document);
    if(average_pf_distance > 0) {
        encode_key_value("average pf distance", average_pf_distance, container, document);
    }
    if(average_pf_confidence > 0) {
        encode_key_value("average pf confidence", average_pf_confidence, container, document);
    }
    encode_key_value("pf fraction", pf_fraction, container, document);
    encode_key_value("pf pooled fraction", pf_pooled_fraction, container, document);
    if(pf_pooled_classified_fraction > 0) {
        encode_key_value("pf pooled classified fraction", pf_pooled_classified_fraction, container, document);
    }
};
AccumulatingTag& AccumulatingTag::operator+=(const AccumulatingTag& rhs) {
    count += rhs.count;
    pf_count += rhs.pf_count;
    accumulated_distance += rhs.accumulated_distance;
    accumulated_confidence += rhs.accumulated_confidence;
    low_conditional_confidence_count += rhs.low_conditional_confidence_count;
    low_confidence_count += rhs.low_confidence_count;
    accumulated_pf_distance += rhs.accumulated_pf_distance;
    accumulated_pf_confidence += rhs.accumulated_pf_confidence;
    return *this;
};

/*  AccumulatingClassifier  */

AccumulatingClassifier::AccumulatingClassifier(const int32_t index) :
    index(index),
    count(0),
    pf_count(0),
    classified_count(0),
    accumulated_classified_distance(0),
    accumulated_classified_confidence(0),
    low_conditional_confidence_count(0),
    low_confidence_count(0),
    pf_classified_count(0),
    accumulated_pf_classified_distance(0),
    accumulated_pf_classified_confidence(0),

    pf_fraction(0),
    classified_fraction(0),
    average_classified_distance(0),
    average_classified_confidence(0),
    pf_classified_fraction(0),
    classified_pf_fraction(0),
    average_pf_classified_distance(0),
    average_pf_classified_confidence(0) {
};
AccumulatingClassifier::AccumulatingClassifier(const AccumulatingClassifier& other) :
    index(other.index),
    count(other.count),
    pf_count(other.pf_count),
    classified_count(other.classified_count),
    accumulated_classified_distance(other.accumulated_classified_distance),
    accumulated_classified_confidence(other.accumulated_classified_confidence),
    low_conditional_confidence_count(other.low_conditional_confidence_count),
    low_confidence_count(other.low_confidence_count),
    pf_classified_count(other.pf_classified_count),
    accumulated_pf_classified_distance(other.accumulated_pf_classified_distance),
    accumulated_pf_classified_confidence(other.accumulated_pf_classified_confidence) {
};
void AccumulatingClassifier::finalize() {
    if(count > 0) {
        pf_fraction = double(pf_count) / double(count);
        classified_fraction = double(classified_count) / double(count);
    }
    if(pf_count > 0) {
        pf_classified_fraction = double(pf_classified_count) / double(pf_count);
    }
    if(classified_count > 0) {
        average_classified_distance = accumulated_classified_distance / double(classified_count);
        average_classified_confidence = accumulated_classified_confidence / double(classified_count);
        classified_pf_fraction = double(pf_classified_count) / double(classified_count);
    }
    if(pf_classified_count > 0) {
        average_pf_classified_distance = accumulated_pf_classified_distance / double(pf_classified_count);
        average_pf_classified_confidence = accumulated_pf_classified_confidence / double(pf_classified_count);
    }
};
void AccumulatingClassifier::encode(Value& container, Document& document) const {
    encode_key_value("index", index, container, document);
    encode_key_value("count", count, container, document);
    encode_key_value("pf count", pf_count, container, document);
    encode_key_value("classified count", classified_count, container, document);

    if(low_conditional_confidence_count > 0) {
        encode_key_value("low conditional confidence count", low_conditional_confidence_count, container, document);
    }
    if(low_confidence_count > 0) {
        encode_key_value("low confidence count", low_confidence_count, container, document);
    }
    encode_key_value("pf classified count", pf_classified_count, container, document);
    encode_key_value("pf fraction", pf_fraction, container, document);
    encode_key_value("classified fraction", classified_fraction, container, document);
    if(average_classified_distance > 0) {
        encode_key_value("average classified distance", average_classified_distance, container, document);
    }
    if(average_classified_confidence > 0) {
        encode_key_value("average classified confidence", average_classified_confidence, container, document);
    }
    encode_key_value("pf classified fraction", pf_classified_fraction, container, document);
    encode_key_value("classified pf fraction", classified_pf_fraction, container, document);
    if(average_pf_classified_distance > 0) {
        encode_key_value("average pf classified distance", average_pf_classified_distance, container, document);
    }
    if(average_pf_classified_confidence > 0) {
        encode_key_value("average pf classified confidence", average_pf_classified_confidence, container, document);
    }
};
AccumulatingClassifier& AccumulatingClassifier::operator+=(const AccumulatingClassifier& rhs) {
    count += rhs.count;
    pf_count += rhs.pf_count;
    classified_count += rhs.classified_count;
    accumulated_classified_distance += rhs.accumulated_classified_distance;
    low_conditional_confidence_count += rhs.low_conditional_confidence_count;
    low_confidence_count += rhs.low_confidence_count;
    accumulated_classified_confidence += rhs.accumulated_classified_confidence;
    pf_classified_count += rhs.pf_classified_count;
    accumulated_pf_classified_distance += rhs.accumulated_pf_classified_distance;
    accumulated_pf_classified_confidence += rhs.accumulated_pf_classified_confidence;
    return *this;
};

/*  NucleotideAccumulator */

NucleotideAccumulator::NucleotideAccumulator() :
    count(0),
    min_quality(0),
    max_quality(0),
    sum_quality(0),
    mean_quality(0),
    Q1(0),
    Q3(0),
    IQR(0),
    LW(0),
    RW(0),
    median_quality(0),
    distribution(EFFECTIVE_PHRED_RANGE, 0) {
};
void NucleotideAccumulator::finalize() {
    for(auto& q : distribution) {
        count += q;
    }
    if(count > 0) {
        for(size_t q(0); q < distribution.size(); ++q) {
            const uint64_t value(distribution[q]);
            sum_quality += (value * q);
            if(value != 0) {
                max_quality = q;
                if(min_quality == 0) {
                    min_quality = q;
                }
            }
        }
        mean_quality = double(sum_quality) / double(count);
        median_quality = quantile(0.5);
        Q1 = quantile(0.25);
        Q3 = quantile(0.75);
        IQR = Q3 - Q1;

        double W(Q1 - IQR * 1.5);
        LW = (W < min_quality) ? min_quality : W;

        W = Q3 + IQR * 1.5;
        RW = (W > max_quality) ? max_quality : W;
    }
};
NucleotideAccumulator& NucleotideAccumulator::operator=(const NucleotideAccumulator& rhs) {
    if(this != &rhs) {
        count = rhs.count;
        min_quality = rhs.min_quality;
        max_quality = rhs.max_quality;
        sum_quality = rhs.sum_quality;
        mean_quality = rhs.mean_quality;
        Q1 = rhs.Q1;
        Q3 = rhs.Q3;
        IQR = rhs.IQR;
        LW = rhs.LW;
        RW = rhs.RW;
        median_quality = rhs.median_quality;
        distribution = rhs.distribution;
    }
    return *this;
};
NucleotideAccumulator& NucleotideAccumulator::operator+=(const NucleotideAccumulator& rhs) {
    for(size_t q(0); q < distribution.size(); ++q) {
        distribution[q] += rhs.distribution[q];
    }
    return *this;
};

/*  CycleAccumulator */

CycleAccumulator::CycleAccumulator() :
    nucleotide_by_code(IUPAC_CODE_SIZE) {
};
void CycleAccumulator::finalize() {
    /* accumulate all nucleotide variations in the NO_NUCLEOTIDE accumulative distribution */
    for(uint8_t i(1); i < nucleotide_by_code.size(); ++i) {
        for(uint8_t p(0); p < EFFECTIVE_PHRED_RANGE; ++p) {
            nucleotide_by_code[NO_NUCLEOTIDE].distribution[p] += nucleotide_by_code[i].distribution[p];
        }
    }
    for(auto& distribution : nucleotide_by_code) {
        distribution.finalize();
    }
};
CycleAccumulator& CycleAccumulator::operator=(const CycleAccumulator& rhs) {
    if(this != &rhs) {
        nucleotide_by_code = rhs.nucleotide_by_code;
    }
    return *this;
};
CycleAccumulator& CycleAccumulator::operator+=(const CycleAccumulator& rhs) {
    for(size_t i(0); i < nucleotide_by_code.size(); ++i) {
        nucleotide_by_code[i] += rhs.nucleotide_by_code[i];
    }
    return *this;
};
