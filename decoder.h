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
#include "transform.h"
#include "channel.h"

class Decoder {
    public:
        Decoder(const Value& ontology) {
        };
        virtual ~Decoder() {
        };
        virtual void decode(const Read& read, Read& output) = 0;
};

template < class T > class DiscreteDecoder : public Decoder {
    public:
        vector< T > element_by_index;
        T undetermined;
        T* decoded;
        DiscreteDecoder(const Value& ontology) :
            Decoder(ontology),
            element_by_index(decode_value_by_key< vector< T > >("codec", ontology)),
            undetermined(find_value_by_key("undetermined", ontology)),
            decoded(NULL) {
        };
};

class ReadGroupDecoder : public DiscreteDecoder< Channel > {
    private:
        string key_buffer;

    protected:
        unordered_map< string, Channel* > channel_by_rg;

    public:
        ReadGroupDecoder(const Value& ontology) :
            DiscreteDecoder< Channel >(ontology) {

            channel_by_rg.reserve(element_by_index.size());
            for(auto& element : element_by_index) {
                channel_by_rg.emplace(make_pair(string(element.rg.ID.s, element.rg.ID.l), &element));
            }
        };
        inline void decode(const Read& input, Read& output) override {
            decoded = &undetermined;
            if(!ks_empty(input.RG())) {
                key_buffer.assign(input.RG().s, input.RG().l);
                auto record = channel_by_rg.find(key_buffer);
                if(record != channel_by_rg.end()) {
                    decoded = record->second;
                }
            }
        };
};

template < class T > class BarcodeDecoder : public DiscreteDecoder< T > {
    protected:
        const Rule rule;
        const int32_t nucleotide_cardinality;
        Observation observation;
        int32_t distance;

    public:
        inline const int32_t segment_cardinality() const {
            return static_cast< int32_t >(observation.segment_cardinality());
        };
        BarcodeDecoder(const Value& ontology) :
            DiscreteDecoder< T >(ontology),
            rule(decode_value_by_key< Rule >("template", ontology)),
            nucleotide_cardinality(decode_value_by_key< int32_t >("nucleotide cardinality", ontology)),
            observation(decode_value_by_key< int32_t >("segment cardinality", ontology)),
            distance(0) {
        };
};

template < class T > class MDDecoder : public BarcodeDecoder< T > {
    protected:
        const uint8_t quality_masking_threshold;
        const vector< uint8_t > distance_tolerance;
        unordered_map< string, T* > element_by_sequence;

    private:
        inline bool match(const T& barcode) {
            bool result(true);
            this->distance = 0;
            if(this->quality_masking_threshold > 0) {
                for(size_t i(0); i < this->observation.segment_cardinality(); ++i) {
                    int32_t error(this->observation[i].masked_distance_from(barcode[i], this->quality_masking_threshold));
                    this->distance += error;
                    if(error > this->distance_tolerance[i]) {
                        result = false;
                    }
                }
            } else {
                for(size_t i(0); i < this->observation.segment_cardinality(); ++i) {
                    int32_t error(this->observation[i].distance_from(barcode[i]));
                    this->distance += error;
                    if(error > this->distance_tolerance[i]) {
                        result = false;
                    }
                }
            }
            return result;
        };

    public:
        MDDecoder(const Value& ontology) :
            BarcodeDecoder< T >(ontology),
            quality_masking_threshold(decode_value_by_key< uint8_t >("quality masking threshold", ontology)),
            distance_tolerance(decode_value_by_key< vector< uint8_t > >("distance tolerance", ontology)) {

            for(auto& element : this->element_by_index) {
                element_by_sequence.emplace(make_pair(string(element), &element));
            }
        };
        inline void decode(const Read& input, Read& output) override {
            this->observation.clear();
            this->decoded = &this->undetermined;
            this->distance = this->nucleotide_cardinality;
            this->rule.apply(input, this->observation);

            /* First try a perfect match to the full barcode sequence */
            auto record = element_by_sequence.find(this->observation);
            if(record != element_by_sequence.end()) {
                this->distance = 0;
                this->decoded = record->second;

            } else {
                /* If no exact match was found try error correction */
                for(auto& barcode : this->element_by_index) {
                    if(match(barcode)) {
                        this->decoded = &barcode;
                        break;
                    }
                }
            }
        };
};

template < class T > class PAMLDecoder : public BarcodeDecoder< T > {
    protected:
        const double noise;
        const double confidence_threshold;
        const double random_barcode_probability;
        const double adjusted_noise_probability;
        double conditioned_decoding_probability;
        double decoding_probability;
        double error_probability;

    public:
        PAMLDecoder(const Value& ontology) :
            BarcodeDecoder< T >(ontology),
            noise(decode_value_by_key< double >("noise", ontology)),
            confidence_threshold(decode_value_by_key< double >("confidence threshold", ontology)),
            random_barcode_probability(1.0 / double(pow(4, (this->nucleotide_cardinality)))),
            adjusted_noise_probability(noise * random_barcode_probability),
            conditioned_decoding_probability(0),
            decoding_probability(0),
            error_probability(1) {
        };
        inline void decode(const Read& input, Read& output) override {
            this->observation.clear();
            this->decoded = &this->undetermined;
            this->distance = this->nucleotide_cardinality;
            conditioned_decoding_probability = 0;
            decoding_probability = 0;
            error_probability = 1;
            this->rule.apply(input, this->observation);

            /*  Compute P(observed|barcode) for each barcode
                Keep track of the channel that yield the maximal prior adjusted probability.
                If r is the observed sequence and b is the barcode sequence
                P(r|b) is the probability that r was observed given b was sequenced.
                Accumulate all prior adjusted probabilities P(b) * P(r|b), in sigma
                using the Kahan summation algorithm to minimize floating point drift
                see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
            */
            double adjusted(0);
            double compensation(0);
            double sigma(0);
            double y(0);
            double t(0);
            double c(0);
            double p(0);
            int32_t d(0);
            for(auto& barcode : this->element_by_index) {
                barcode.accurate_decoding_probability(this->observation, c, d);
                p = c * barcode.concentration;
                y = p - compensation;
                t = sigma + y;
                compensation = (t - sigma) - y;
                sigma = t;
                if(p > adjusted) {
                    this->decoded = &barcode;
                    conditioned_decoding_probability = c;
                    this->distance = d;
                    adjusted = p;
                }
            }
            /*  Compute P(barcode|observed)
                P(b|r), the probability that b was sequenced given r was observed
                P(b|r) = P(r|b) * P(b) / ( P(noise) * P(r|noise) + sigma )
                where sigma = sum of P(r|b) * P(b) over b */
            decoding_probability = adjusted / (sigma + adjusted_noise_probability);

            /* Check for decoding failure and assign to the undetermined channel if decoding failed */
            if(conditioned_decoding_probability > random_barcode_probability && decoding_probability > confidence_threshold) {
                error_probability = 1 - decoding_probability;
            } else {
                this->decoded = &this->undetermined;
            }
        };
};

class MultiplexMDDecoder : public MDDecoder< Channel > {
    public:
        MultiplexMDDecoder(const Value& ontology) :
            MDDecoder< Channel >(ontology) {
        };
        inline void decode(const Read& input, Read& output) override {
            MDDecoder< Channel >::decode(input, output);
            output.assign_RG(this->decoded->rg);
            output.set_multiplex_barcode(this->observation);
            output.set_multiplex_distance(this->distance);
        };
};

class MultiplexPAMLDecoder : public PAMLDecoder< Channel > {
    public:
        MultiplexPAMLDecoder(const Value& ontology) :
            PAMLDecoder< Channel >(ontology) {
        };
        inline void decode(const Read& input, Read& output) override {
            PAMLDecoder< Channel >::decode(input, output);
            output.assign_RG(this->decoded->rg);
            output.set_multiplex_barcode(this->observation);
            output.set_multiplex_error_probability(this->error_probability);
            output.set_multiplex_distance(this->distance);
        };
};

class CellularMDDecoder : public MDDecoder< Barcode > {
    public:
        CellularMDDecoder(const Value& ontology) :
            MDDecoder< Barcode >(ontology) {
        };
        inline void decode(const Read& input, Read& output) override {
            MDDecoder< Barcode >::decode(input, output);
            output.update_raw_cellular_barcode(this->observation);
            output.update_cellular_barcode(*this->decoded);
            if(this->decoded != &this->undetermined) {
                output.update_cellular_distance(this->distance);
            } else {
                output.set_cellular_distance(0);
            }
        };
};

class CellularPAMLDecoder : public PAMLDecoder< Barcode > {
    public:
        CellularPAMLDecoder(const Value& ontology) :
            PAMLDecoder< Barcode >(ontology) {
        };
        inline void decode(const Read& input, Read& output) override {
            PAMLDecoder< Barcode >::decode(input, output);
            output.update_raw_cellular_barcode(this->observation);
            output.update_cellular_barcode(*this->decoded);
            if(this->decoded != &this->undetermined) {
                output.update_cellular_error_probability(this->error_probability);
                output.update_cellular_distance(this->distance);
            } else {
                output.set_cellular_error_probability(1);
                output.set_cellular_distance(0);
            }
        };
};

class MolecularSimpleDecoder : public Decoder {
    protected:
        const int32_t nucleotide_cardinality;
        const Rule rule;
        Observation observation;

    public:
        MolecularSimpleDecoder(const Value& ontology) :
            Decoder(ontology),
            nucleotide_cardinality(decode_value_by_key< int32_t >("nucleotide cardinality", ontology)),
            rule(decode_value_by_key< Rule >("template", ontology)),
            observation(decode_value_by_key< int32_t >("segment cardinality", ontology)) {
        };
        inline void decode(const Read& input, Read& output) override {
            observation.clear();
            rule.apply(input, observation);
            output.update_molecular_barcode(observation);
        };
};

#endif /* PHENIQS_DECODER_H */
