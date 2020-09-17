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

#include "selector.h"

/*  AccumulatingOption */

AccumulatingOption::AccumulatingOption() :
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
    pf_pooled_classified_fraction(0),
    estimated_concentration_prior(0) {
};
AccumulatingOption::AccumulatingOption(const AccumulatingOption& other) :
    count(other.count),
    pf_count(other.pf_count),
    accumulated_distance(other.accumulated_distance),
    accumulated_confidence(other.accumulated_confidence),
    low_conditional_confidence_count(other.low_conditional_confidence_count),
    low_confidence_count(other.low_confidence_count),
    accumulated_pf_distance(other.accumulated_pf_distance),
    accumulated_pf_confidence(other.accumulated_pf_confidence),

    pf_fraction(0),
    average_distance(0),
    average_confidence(0),
    average_pf_distance(0),
    average_pf_confidence(0),
    pooled_fraction(0),
    pf_pooled_fraction(0),
    pooled_classified_fraction(0),
    pf_pooled_classified_fraction(0),
    estimated_concentration_prior(0) {
};
void AccumulatingOption::collect(const AccumulatingOption& other) {
    count += other.count;
    pf_count += other.pf_count;
    accumulated_distance += other.accumulated_distance;
    accumulated_confidence += other.accumulated_confidence;
    low_conditional_confidence_count += other.low_conditional_confidence_count;
    low_confidence_count += other.low_confidence_count;
    accumulated_pf_distance += other.accumulated_pf_distance;
    accumulated_pf_confidence += other.accumulated_pf_confidence;
};
void AccumulatingOption::finalize(const AccumulatingSelector& parent) {
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
void AccumulatingOption::encode(Value& container, Document& document) const {
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
    if(estimated_concentration_prior > 0) {
        encode_key_value("estimated concentration prior", estimated_concentration_prior, container, document);
    }
};

/*  AccumulatingSelector  */

AccumulatingSelector::AccumulatingSelector(const int32_t index) :
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
    average_pf_classified_confidence(0),
    estimated_noise_prior(0) {
};
AccumulatingSelector::AccumulatingSelector(const AccumulatingSelector& other) :
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
    accumulated_pf_classified_confidence(other.accumulated_pf_classified_confidence),

    pf_fraction(0),
    classified_fraction(0),
    average_classified_distance(0),
    average_classified_confidence(0),
    pf_classified_fraction(0),
    classified_pf_fraction(0),
    average_pf_classified_distance(0),
    average_pf_classified_confidence(0),
    estimated_noise_prior(0) {
};
void AccumulatingSelector::collect(const AccumulatingSelector& other) {
    count += other.count;
    pf_count += other.pf_count;
    classified_count += other.classified_count;
    accumulated_classified_distance += other.accumulated_classified_distance;
    low_conditional_confidence_count += other.low_conditional_confidence_count;
    low_confidence_count += other.low_confidence_count;
    accumulated_classified_confidence += other.accumulated_classified_confidence;
    pf_classified_count += other.pf_classified_count;
    accumulated_pf_classified_distance += other.accumulated_pf_classified_distance;
    accumulated_pf_classified_confidence += other.accumulated_pf_classified_confidence;
};
void AccumulatingSelector::finalize() {
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
void AccumulatingSelector::encode(Value& container, Document& document) const {
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
    if(estimated_noise_prior > 0) {
        encode_key_value("estimated noise prior", estimated_noise_prior, container, document);
    }
};
