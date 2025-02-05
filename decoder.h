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

#ifndef PHENIQS_DECODER_H
#define PHENIQS_DECODER_H

#include "include.h"
#include "classifier.h"
#include "transform.h"

template < class T > class Decoder : public Classifier< T > {
    protected:
        const Rule rule;
        const int32_t nucleotide_cardinality;
        const uint8_t high_quality_threshold;
        const int32_t high_quality_distance_threshold;
        Observation observation;
        int32_t edit_distance;
        int32_t high_quality_edit_distance;


    public:
        inline const int32_t segment_cardinality() const {
            return static_cast< int32_t >(observation.segment_cardinality());
        };
        Decoder(const Value& ontology) try :
            Classifier< T >(ontology),
            rule(decode_value_by_key< Rule >("transform", ontology)),
            nucleotide_cardinality(decode_value_by_key< int32_t >("nucleotide cardinality", ontology)),
            high_quality_threshold(decode_value_by_key< uint8_t >("high quality threshold", ontology)),
            high_quality_distance_threshold(decode_value_by_key< int32_t >("high quality distance threshold", ontology)),
            observation(decode_value_by_key< int32_t >("segment cardinality", ontology)),
            edit_distance(0),
            high_quality_edit_distance(0) {

            } catch(Error& error) {
                error.push("Decoder");
                throw;
        };
        Decoder(const Decoder< T >& other) :
            Classifier< T >(other),
            rule(other.rule),
            nucleotide_cardinality(other.nucleotide_cardinality),
            high_quality_threshold(other.high_quality_threshold),
            high_quality_distance_threshold(other.high_quality_distance_threshold),
            observation(other.observation.segment_cardinality()),
            edit_distance(0),
            high_quality_edit_distance(0) {
        };
        inline void classify(const Read& input, Read& output) override {
            if(this->decoded->is_classified() && edit_distance) {
                this->decoded->accumulated_distance += static_cast< uint64_t >(edit_distance);
                if(!output.qcfail()) {
                    this->decoded->accumulated_pf_distance += static_cast< uint64_t >(edit_distance);
                }
            }
            Classifier< T >::classify(input, output);
        };
        inline void finalize() override {
            for(auto& element : this->tag_array) {
                this->accumulated_classified_distance += element.accumulated_distance;
                this->accumulated_pf_classified_distance += element.accumulated_pf_distance;
            }
            Classifier< T >::finalize();
        };
};

#endif /* PHENIQS_DECODER_H */
