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

#ifndef PHENIQS_BARCODE_H
#define PHENIQS_BARCODE_H

#include "include.h"
#include "sequence.h"

class Barcode : public SequenceArray< Sequence > {
    friend ostream& operator<<(ostream& o, const Barcode& barcode);
    friend bool encode_key_value(const string& key, const Barcode& value, Value& node, Document& document);
    void operator=(Barcode const &) = delete;

    public:
        const int32_t index;
        const double concentration;
        Barcode(const Value& ontology);
        Barcode(const Barcode& other) :
            SequenceArray< Sequence >(other),
            index(other.index),
            concentration(other.concentration) {
        };
        operator string() const {
            /* NOTICE this is in BAM encoding not iupac and will not look as expected when printed */
            string key;
            for(const auto& segment : segment_array) {
                for(int32_t i(0); i < segment.length; ++i) {
                    key.push_back(segment.code[i]);
                }
            }
            return key;
        };
        inline string iupac_ambiguity() const {
            string value;
            for(const auto& segment : segment_array) {
                for(int32_t i = 0; i < segment.length; ++i) {
                    value.push_back(BamToAmbiguousAscii[segment.code[i]]);
                }
            }
            return value;
        };
        inline void decoding_probability(const Observation& observation, double& probability, int32_t& distance) const {
            double p(1);
            int32_t d(0);
            for(size_t i(0); i < segment_array.size(); ++i) {
                const Sequence& reference = segment_array[i];
                const ObservedSequence& observed = observation[i];
                for(int32_t j = 0; j < reference.length; ++j) {
                    if(observed.code[j] == reference.code[j]) {
                        p *= quality_to_inverse_probability(observed.quality[j]);
                    } else {
                        d += 1;
                        if(observed.code[j] != ANY_NUCLEOTIDE) {
                            p *= quality_to_third_probability(observed.quality[j]);
                        } else {
                            p *= UNIFORM_BASE_PROBABILITY;
                        }
                    }
                }
            }
            distance = d;
            probability = p;
        };
        inline void accurate_decoding_probability(const Observation& observation, double& probability, int32_t& distance) const {
            double q(0);
            int32_t d(0);
            for(size_t i(0); i < segment_array.size(); ++i) {
                const Sequence& reference = segment_array[i];
                const ObservedSequence& observed = observation[i];
                for(int32_t j(0); j < reference.length; ++j) {
                    if(observed.code[j] == reference.code[j]) {
                        q += quality_to_inverse_quality(observed.quality[j]);
                    } else {
                        d += 1;
                        if(observed.code[j] != ANY_NUCLEOTIDE) {
                            q += double(observed.quality[j]);
                        } else {
                            q += UNIFORM_BASE_PHRED;
                        }
                    }
                }
            }
            distance = d;
            probability = pow(10.0, q * -0.1);
        };
        inline void compensated_decoding_probability(const Observation& observation, double& probability, int32_t& distance) const {
            // use the Kahan summation algorithm to minimize floating point drift
            // see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
            double sigma(0);
            double compensation(0);
            double y(0);
            double t(0);
            double q(0);
            int32_t d(0);
            for(size_t i(0); i < segment_array.size(); ++i) {
                const Sequence& reference = segment_array[i];
                const ObservedSequence& observed = observation[i];
                for(int32_t j(0); j < reference.length; ++j) {
                    if(observed.code[j] == reference.code[j]) {
                        q = quality_to_inverse_quality(observed.quality[j]);
                    } else {
                        d += 1;
                        if(observed.code[j] != ANY_NUCLEOTIDE) {
                            q = double(observed.quality[j]);
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
};

ostream& operator<<(ostream& o, const Barcode& barcode);
template<> vector< Barcode > decode_value_by_key(const Value::Ch* key, const Value& container);
bool encode_key_value(const string& key, const Barcode& value, Value& container, Document& document);

#endif /* PHENIQS_BARCODE_H */
