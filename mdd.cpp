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

#include "mdd.h"

template < class T > MinimumDistanceDecoder< T >::MinimumDistanceDecoder(const Value& ontology) try :
    ObservingDecoder< T >(ontology),
    quality_masking_threshold(decode_value_by_key< uint8_t >("quality masking threshold", ontology)),
    distance_tolerance(decode_value_by_key< vector< int32_t > >("distance tolerance", ontology)) {

    for(auto& element : this->element_by_index) {
        element_by_sequence.emplace(make_pair(string(element), &element));
    }

    } catch(Error& error) {
        error.push("MinimumDistanceDecoder");
        throw;
};
template < class T > void MinimumDistanceDecoder< T >::decode(const Read& input, Read& output) {
    this->observation.clear();
    this->rule.apply(input, this->observation);
    this->decoded = &this->unclassified;
    this->decoding_hamming_distance = 0;

    /* First try a perfect match to the full barcode sequence */
    auto record = element_by_sequence.find(this->observation);
    if(record != element_by_sequence.end()) {
        this->decoded = record->second;

    } else {
        /* If no exact match was not found try error correction */
        for(auto& barcode : this->element_by_index) {
            int32_t distance(0);
            bool successful(true);
            if(this->quality_masking_threshold > 0) {
                for(size_t i(0); i < this->observation.segment_cardinality(); ++i) {
                    int32_t error(this->observation[i].masked_distance_from(barcode[i], this->quality_masking_threshold));
                    if(error > this->distance_tolerance[i]) {
                        successful = false;
                        break;
                    } else {
                        distance += error;
                    }
                }
            } else {
                for(size_t i(0); i < this->observation.segment_cardinality(); ++i) {
                    int32_t error(this->observation[i].distance_from(barcode[i]));
                    if(error > this->distance_tolerance[i]) {
                        successful = false;
                        break;
                    } else {
                        distance += error;
                    }
                }
            }

            if(successful) {
                this->decoding_hamming_distance = distance;
                this->decoded = &barcode;
                break;
            }
        }
    }
    ObservingDecoder< T >::decode(input, output);
};

MDMultiplexDecoder::MDMultiplexDecoder(const Value& ontology) try :
    MinimumDistanceDecoder< Channel >(ontology) {

    } catch(Error& error) {
        error.push("MDMultiplexDecoder");
        throw;
};
void MDMultiplexDecoder::decode(const Read& input, Read& output) {
    MinimumDistanceDecoder< Channel >::decode(input, output);
    output.assign_RG(this->decoded->rg);
    output.update_multiplex_barcode(this->observation);
    output.update_multiplex_distance(this->decoding_hamming_distance);
};

MDCellularDecoder::MDCellularDecoder(const Value& ontology) try :
    MinimumDistanceDecoder< Barcode >(ontology) {

    } catch(Error& error) {
        error.push("MDCellularDecoder");
        throw;
};
void MDCellularDecoder::decode(const Read& input, Read& output) {
    MinimumDistanceDecoder< Barcode >::decode(input, output);
    output.update_raw_cellular_barcode(this->observation);
    output.update_cellular_barcode(*this->decoded);
    if(this->decoded->is_classified()) {
        output.update_cellular_distance(this->decoding_hamming_distance);
    } else {
        output.set_cellular_distance(0);
    }
};
