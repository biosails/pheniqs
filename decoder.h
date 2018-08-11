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
        Decoder(const Value& ontology) {};
        virtual ~Decoder() {};
        virtual void decode(const Read& read, Read& output) {};
};

template < class T > class RoutingDecoder : public Decoder {
    public:
        T unclassified;
        T* decoded;
        RoutingDecoder(const Value& ontology) try :
            Decoder(ontology),
            unclassified(find_value_by_key("undetermined", ontology)),
            decoded(NULL) {

            decoded = &unclassified;

            } catch(Error& error) {
                error.push("RoutingDecoder");
                throw;
        };
        void decode(const Read& read, Read& output) override {};

};

template < class T > class TransparentDecoder : public RoutingDecoder< T > {
    public:
        TransparentDecoder(const Value& ontology) try :
            RoutingDecoder< T >(ontology) {

            } catch(Error& error) {
                error.push("TransparentDecoder");
                throw;
        };
};

template < class T > class DiscreteDecoder : public RoutingDecoder< T > {
    public:
        vector< T > element_by_index;
        DiscreteDecoder(const Value& ontology) try :
            RoutingDecoder< T >(ontology),
            element_by_index(decode_value_by_key< vector< T > >("codec", ontology)) {

            } catch(Error& error) {
                error.push("DiscreteDecoder");
                throw;
        };
};

template < class T > class ReadGroupDecoder : public DiscreteDecoder< T > {
    public:
        ReadGroupDecoder(const Value& ontology) try :
            DiscreteDecoder< T >(ontology) {

            element_by_rg.reserve(this->element_by_index.size());
            for(auto& element : this->element_by_index) {
                element_by_rg.emplace(make_pair(string(element.rg.ID.s, element.rg.ID.l), &element));
            }

            } catch(Error& error) {
                error.push("ReadGroupDecoder");
                throw;
        };
        inline void decode(const Read& input, Read& output) override {
            this->decoded = &this->unclassified;
            if(!ks_empty(input.RG())) {
                rg_id_buffer.assign(input.RG().s, input.RG().l);
                auto record = element_by_rg.find(rg_id_buffer);
                if(record != element_by_rg.end()) {
                    this->decoded = record->second;
                }
            }
        };

    protected:
        unordered_map< string, T* > element_by_rg;

    private:
        string rg_id_buffer;

};

template < class T > class ObservationDecoder : public DiscreteDecoder< T > {
    protected:
        const Rule rule;
        const int32_t nucleotide_cardinality;
        Observation observation;
        int32_t decoding_distance;

    public:
        inline const int32_t segment_cardinality() const {
            return static_cast< int32_t >(observation.segment_cardinality());
        };
        ObservationDecoder(const Value& ontology) try :
            DiscreteDecoder< T >(ontology),
            rule(decode_value_by_key< Rule >("transform", ontology)),
            nucleotide_cardinality(decode_value_by_key< int32_t >("nucleotide cardinality", ontology)),
            observation(decode_value_by_key< int32_t >("segment cardinality", ontology)),
            decoding_distance(0) {

            } catch(Error& error) {
                error.push("ObservationDecoder");
                throw;
        };
};

template < class T > class MDDecoder : public ObservationDecoder< T > {
    protected:
        const uint8_t quality_masking_threshold;
        const vector< int32_t > distance_tolerance;
        unordered_map< string, T* > element_by_sequence;

    public:
        MDDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;

    private:
        inline bool match(T& barcode);
};

template < class T > class PAMLDecoder : public ObservationDecoder< T > {
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

class MolecularNaiveDecoder : public Decoder {
    protected:
        const int32_t nucleotide_cardinality;
        const Rule rule;
        Observation observation;

    public:
        MolecularNaiveDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

#endif /* PHENIQS_DECODER_H */
