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
#include "accumulator.h"

class Barcode : public SequenceArray< Sequence >, public AccumulatingIdentifier {
    friend ostream& operator<<(ostream& o, const Barcode& barcode);
    friend bool encode_key_value(const string& key, const Barcode& value, Value& node, Document& document);
    void operator=(Barcode const &) = delete;

    public:
        const int32_t index;
        const double concentration;
        Barcode(const Value& ontology);
        Barcode(const Barcode& other);
        inline const bool is_classified() const {
            /* by convention, enforced by the job configuration loader, barcode 0 is always the unclassified */
            return index > 0;
        };
        inline const bool is_unclassified() const {
            return index == 0;
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
        inline void decoding_probability(const Observation& observation, double& probability, int32_t& distance) const {
            probability = 1;
            for(size_t i(0); i < segment_array.size(); ++i) {
                const Sequence& sequnced = segment_array[i];
                const ObservedSequence& observed = observation[i];
                for(int32_t j(0); j < sequnced.length; ++j) {
                    if(observed.code[j] == sequnced.code[j]) {
                        probability *= quality_to_inverse_probability(observed.quality[j]);
                    } else {
                        ++distance;
                        if(observed.code[j] == ANY_NUCLEOTIDE) {
                            probability *= UNIFORM_BASE_PROBABILITY;
                        } else { probability *= quality_to_third_probability(observed.quality[j]); }
                    }
                }
            }
        };
        inline void accurate_decoding_probability(const Observation& observation, double& probability) const {
            double sigma_q(0);
            uint32_t sigma_c(0);
            for(size_t i(0); i < segment_array.size(); ++i) {
                const Sequence& sequnced = segment_array[i];
                const ObservedSequence& observed = observation[i];
                for(int32_t j(0); j < sequnced.length; ++j) {
                    if(observed.code[j] == sequnced.code[j]) {
                        sigma_q += quality_to_inverse_quality(observed.quality[j]);
                    } else {
                        if(observed.code[j] == ANY_NUCLEOTIDE) {
                            sigma_q += UNIFORM_BASE_PHRED;
                        } else { sigma_c += observed.quality[j]; }
                    }
                }
            }
            sigma_q += sigma_c;
            probability = pow(PHRED_PROBABILITY_FACTOR, sigma_q);
        };
        inline void accurate_decoding_probability(const Observation& observation, double& probability, int32_t& distance) const {
            /*  sigma_c accumulates uint8_t Phred quality scores
                sigma_q accumulates double precision floats from inversed quality scores and UNIFORM_BASE_PHRED
            */
            distance = 0;
            double sigma_q(0);
            uint32_t sigma_c(0);
            for(size_t i(0); i < segment_array.size(); ++i) {
                const Sequence& sequnced = segment_array[i];
                const ObservedSequence& observed = observation[i];
                for(int32_t j(0); j < sequnced.length; ++j) {
                    if(observed.code[j] == sequnced.code[j]) {
                        sigma_q += quality_to_inverse_quality(observed.quality[j]);
                    } else {
                        ++distance;
                        if(observed.code[j] == ANY_NUCLEOTIDE) {
                            sigma_q += UNIFORM_BASE_PHRED;
                        } else { sigma_c += observed.quality[j]; }
                    }
                }
            }
            sigma_q += sigma_c;
            probability = pow(PHRED_PROBABILITY_FACTOR, sigma_q);
        };
        inline void compensated_decoding_probability(const Observation& observation, double& probability) const {
            /*  use the Kahan summation algorithm to minimize floating point drift
                see https://en.wikipedia.org/wiki/Kahan_summation_algorithm

                sigma_c accumulates uint8_t Phred quality scores
                sigma_q accumulates double precision floats from inversed quality scores and UNIFORM_BASE_PHRED
            */
            double y(0);
            double t(0);
            double sigma_q(0);
            uint32_t sigma_c(0);
            double compensation(0);
            for(size_t i(0); i < segment_array.size(); ++i) {
                const Sequence& sequnced = segment_array[i];
                const ObservedSequence& observed = observation[i];
                for(int32_t j(0); j < sequnced.length; ++j) {
                    if(observed.code[j] == sequnced.code[j]) {
                        y = quality_to_inverse_quality(observed.quality[j]) - compensation;
                        t = sigma_q + y;
                        compensation = (t - sigma_q) - y;
                        sigma_q = t;

                    } else {
                        if(observed.code[j] == ANY_NUCLEOTIDE) {
                            y = UNIFORM_BASE_PHRED - compensation;
                            t = sigma_q + y;
                            compensation = (t - sigma_q) - y;
                            sigma_q = t;
                        } else { sigma_c += observed.quality[j]; }
                    }
                }
            }
            sigma_q += sigma_c;
            probability = pow(PHRED_PROBABILITY_FACTOR, sigma_q);
        };
        inline void compensated_decoding_probability(const Observation& observation, double& probability, int32_t& distance) const {
            /*  use the Kahan summation algorithm to minimize floating point drift
                see https://en.wikipedia.org/wiki/Kahan_summation_algorithm

                sigma_c accumulates uint8_t Phred quality scores
                sigma_q accumulates double precision floats from inversed quality scores and UNIFORM_BASE_PHRED
            */
            double y(0);
            double t(0);
            distance = 0;
            double sigma_q(0);
            uint32_t sigma_c(0);
            double compensation(0);
            for(size_t i(0); i < segment_array.size(); ++i) {
                const Sequence& sequnced = segment_array[i];
                const ObservedSequence& observed = observation[i];
                for(int32_t j(0); j < sequnced.length; ++j) {
                    if(observed.code[j] == sequnced.code[j]) {
                        y = quality_to_inverse_quality(observed.quality[j]) - compensation;
                        t = sigma_q + y;
                        compensation = (t - sigma_q) - y;
                        sigma_q = t;

                    } else {
                        ++distance;
                        if(observed.code[j] == ANY_NUCLEOTIDE) {
                            y = UNIFORM_BASE_PHRED - compensation;
                            t = sigma_q + y;
                            compensation = (t - sigma_q) - y;
                            sigma_q = t;
                        } else { sigma_c += observed.quality[j]; }
                    }
                }
            }
            sigma_q += sigma_c;
            probability = pow(PHRED_PROBABILITY_FACTOR, sigma_q);
        };
        void encode(Value& container, Document& document) const override;
        Barcode& operator+=(const Barcode& rhs);
};

ostream& operator<<(ostream& o, const Barcode& barcode);
template<> vector< Barcode > decode_value_by_key(const Value::Ch* key, const Value& container);

#endif /* PHENIQS_BARCODE_H */
