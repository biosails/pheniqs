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
#include "accumulate.h"

template < class T > class RoutingDecoder : public AccumulatingDecoder {
    public:
        T* decoded;
        T unclassified;
        vector< T > element_by_index;

        RoutingDecoder(const Value& ontology) try :
            AccumulatingDecoder(),
            decoded(NULL),
            unclassified(find_value_by_key("undetermined", ontology)),
            element_by_index(decode_value_by_key< vector< T > >("codec", ontology)) {

            decoded = &unclassified;

            } catch(Error& error) {
                error.push("RoutingDecoder");
                throw;
        };
        virtual inline void decode(const Read& input, Read& output) {
            ++(decoded->count);
            if(!input.qcfail()) {
                ++(decoded->pf_count);
            }
        };
        inline void finalize() override {
            for(auto& element : element_by_index) {
                this->classified_count += element.count;
                this->pf_classified_count += element.pf_count;
            }
            this->count += this->classified_count + unclassified.count;
            this->pf_count += this->pf_classified_count + unclassified.pf_count;

            for(auto& element : element_by_index) {
                element.finalize(*this);
            }
            unclassified.finalize(*this);
            AccumulatingDecoder::finalize();
        };
        RoutingDecoder< T >& operator+=(const RoutingDecoder< T >& rhs) {
            AccumulatingDecoder::operator+=(rhs);
            unclassified += rhs.unclassified;
            for(size_t index(0); index < element_by_index.size(); ++index) {
                element_by_index[index] += rhs.element_by_index[index];
            }
            return *this;
        };
        void encode(Value& container, Document& document) const override {
            AccumulatingDecoder::encode(container, document);

            Value unclassified_report(kObjectType);
            unclassified.encode(unclassified_report, document);
            container.AddMember("unclassified", unclassified_report.Move(), document.GetAllocator());

            if(!element_by_index.empty()) {
                Value element_report_array(kArrayType);
                for(auto& element : element_by_index) {
                    Value element_report(kObjectType);
                    element.encode(element_report, document);
                    element_report_array.PushBack(element_report.Move(), document.GetAllocator());
                }
                container.AddMember("classified", element_report_array.Move(), document.GetAllocator());
            }
        };
};

template < class T > class ReadGroupDecoder : public RoutingDecoder< T > {
    public:
        ReadGroupDecoder(const Value& ontology) try :
            RoutingDecoder< T >(ontology) {

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
            if(ks_not_empty(input.RG())) {
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

template < class T > class ObservingDecoder : public RoutingDecoder< T > {
    protected:
        const Rule rule;
        const int32_t nucleotide_cardinality;
        Observation observation;
        int32_t decoding_hamming_distance;

    public:
        inline const int32_t segment_cardinality() const {
            return static_cast< int32_t >(observation.segment_cardinality());
        };
        ObservingDecoder(const Value& ontology) try :
            RoutingDecoder< T >(ontology),
            rule(decode_value_by_key< Rule >("transform", ontology)),
            nucleotide_cardinality(decode_value_by_key< int32_t >("nucleotide cardinality", ontology)),
            observation(decode_value_by_key< int32_t >("segment cardinality", ontology)),
            decoding_hamming_distance(0) {

            } catch(Error& error) {
                error.push("ObservingDecoder");
                throw;
        };
        inline void decode(const Read& input, Read& output) override {
            if(this->decoded->is_classified() && decoding_hamming_distance) {
                this->decoded->accumulated_distance += static_cast< uint64_t >(decoding_hamming_distance);
                if(!input.qcfail()) {
                    this->decoded->accumulated_pf_distance += static_cast< uint64_t >(decoding_hamming_distance);
                }
            }
            RoutingDecoder< T >::decode(input, output);
        };
        inline void finalize() override {
            for(auto& element : this->element_by_index) {
                this->accumulated_classified_distance += element.accumulated_distance;
                this->accumulated_pf_classified_distance += element.accumulated_pf_distance;
            }
            RoutingDecoder< T >::finalize();
        };
};

template < class T > class MinimumDistanceDecoder : public ObservingDecoder< T > {
    protected:
        const uint8_t quality_masking_threshold;
        const vector< int32_t > distance_tolerance;
        unordered_map< string, T* > element_by_sequence;

    public:
        MinimumDistanceDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;

    private:
        inline bool match(T& barcode);
};

template < class T > class PhredAdjustedMaximumLikelihoodDecoder : public ObservingDecoder< T > {
    protected:
        const double noise;
        const double confidence_threshold;
        const double random_barcode_probability;
        const double adjusted_noise_probability;
        double conditional_decoding_probability;
        double decoding_confidence;

    public:
        PhredAdjustedMaximumLikelihoodDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
        inline void finalize() override {
            for(auto& element : this->element_by_index) {
                this->accumulated_classified_confidence += element.accumulated_confidence;
                this->accumulated_pf_classified_confidence += element.accumulated_pf_confidence;
            }
            ObservingDecoder< T >::finalize();
        };
};

/* Multiplex */
class MDMultiplexDecoder : public MinimumDistanceDecoder< Channel > {
    public:
        MDMultiplexDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

class PAMLMultiplexDecoder : public PhredAdjustedMaximumLikelihoodDecoder< Channel > {
    public:
        PAMLMultiplexDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

/* Molecular */
class NaiveMolecularDecoder : public ObservingDecoder< Barcode > {
    public:
        NaiveMolecularDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

/* Cellular */
class MDCellularDecoder : public MinimumDistanceDecoder< Barcode > {
    public:
        MDCellularDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

class PAMLCellularDecoder : public PhredAdjustedMaximumLikelihoodDecoder< Barcode > {
    public:
        PAMLCellularDecoder(const Value& ontology);
        inline void decode(const Read& input, Read& output) override;
};

#endif /* PHENIQS_DECODER_H */
