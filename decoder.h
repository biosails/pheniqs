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
        ReadGroupDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

template < class T > class BarcodeDecoder : public DiscreteDecoder< T > {
    protected:
        const Rule rule;
        const int32_t nucleotide_cardinality;
        Observation observation;
        int32_t decoding_distance;

    public:
        inline const int32_t segment_cardinality() const {
            return static_cast< int32_t >(observation.segment_cardinality());
        };
        BarcodeDecoder(const Value& ontology) :
            DiscreteDecoder< T >(ontology),
            rule(decode_value_by_key< Rule >("template", ontology)),
            nucleotide_cardinality(decode_value_by_key< int32_t >("nucleotide cardinality", ontology)),
            observation(decode_value_by_key< int32_t >("segment cardinality", ontology)),
            decoding_distance(0) {
        };
};

template < class T > class MDDecoder : public BarcodeDecoder< T > {
    protected:
        const uint8_t quality_masking_threshold;
        const vector< uint8_t > distance_tolerance;
        unordered_map< string, T* > element_by_sequence;

    public:
        MDDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;

    private:
        inline bool match(T& barcode);
};

template < class T > class PAMLDecoder : public BarcodeDecoder< T > {
    protected:
        const double noise;
        const double confidence_threshold;
        const double random_barcode_probability;
        const double adjusted_noise_probability;
        double conditioned_decoding_probability;
        double decoding_probability;

    public:
        PAMLDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

class MultiplexMDDecoder : public MDDecoder< Channel > {
    public:
        MultiplexMDDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

class MultiplexPAMLDecoder : public PAMLDecoder< Channel > {
    public:
        MultiplexPAMLDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

class CellularMDDecoder : public MDDecoder< Barcode > {
    public:
        CellularMDDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

class CellularPAMLDecoder : public PAMLDecoder< Barcode > {
    public:
        CellularPAMLDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

class MolecularSimpleDecoder : public Decoder {
    protected:
        const int32_t nucleotide_cardinality;
        const Rule rule;
        Observation observation;

    public:
        MolecularSimpleDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

#endif /* PHENIQS_DECODER_H */
