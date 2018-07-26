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

#include "decoder.h"

template < class T > MDDecoder< T >::MDDecoder(const Value& ontology) try :
    ObservationDecoder< T >(ontology),
    quality_masking_threshold(decode_value_by_key< uint8_t >("quality masking threshold", ontology)),
    distance_tolerance(decode_value_by_key< vector< int32_t > >("distance tolerance", ontology)) {

    for(auto& element : this->element_by_index) {
        element_by_sequence.emplace(make_pair(string(element), &element));
    }

    } catch(ConfigurationError& error) {
        throw ConfigurationError("MDDecoder :: " + error.message);

    } catch(exception& error) {
        throw InternalError("MDDecoder :: " + string(error.what()));
};
template < class T > bool MDDecoder< T >::match(T& barcode) {
    bool result(true);
    int32_t hamming_distance(0);

    if(this->quality_masking_threshold > 0) {
        for(size_t i(0); i < this->observation.segment_cardinality(); ++i) {
            int32_t error(this->observation[i].masked_distance_from(barcode[i], this->quality_masking_threshold));
            if(error > this->distance_tolerance[i]) {
                result = false;
                break;
            } else {
                hamming_distance += error;
            }
        }
    } else {
        for(size_t i(0); i < this->observation.segment_cardinality(); ++i) {
            int32_t error(this->observation[i].distance_from(barcode[i]));
            if(error > this->distance_tolerance[i]) {
                result = false;
                break;
            } else {
                hamming_distance += error;
            }
        }
    }

    if(result) {
        this->decoding_distance = hamming_distance;
        this->decoded = &barcode;
    }
    return result;
};
template < class T > void MDDecoder< T >::decode(const Read& input, Read& output) {
    this->observation.clear();
    this->decoded = &this->unclassified;
    this->decoding_distance = 0;
    this->rule.apply(input, this->observation);

    /* First try a perfect match to the full barcode sequence */
    auto record = element_by_sequence.find(this->observation);
    if(record != element_by_sequence.end()) {
        this->decoding_distance = 0;
        this->decoded = record->second;

    } else {
        /* If no exact match was not found try error correction */
        for(auto& barcode : this->element_by_index) {
            if(match(barcode)) {
                break;
            }
        }
    }
};

template < class T > PAMLDecoder< T >::PAMLDecoder(const Value& ontology) try :
    ObservationDecoder< T >(ontology),
    noise(decode_value_by_key< double >("noise", ontology)),
    confidence_threshold(decode_value_by_key< double >("confidence threshold", ontology)),
    random_barcode_probability(1.0 / double(pow(4, (this->nucleotide_cardinality)))),
    adjusted_noise_probability(noise * random_barcode_probability),
    conditioned_decoding_probability(0),
    decoding_probability(0) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("PAMLDecoder :: " + error.message);

    } catch(exception& error) {
        throw InternalError("PAMLDecoder :: " + string(error.what()));
};
template < class T > void PAMLDecoder< T >::decode(const Read& input, Read& output) {
    this->observation.clear();
    this->decoded = &this->unclassified;
    this->decoding_distance = 0;
    this->decoding_probability = 0;
    this->conditioned_decoding_probability = 0;
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
            this->decoding_distance = d;
            adjusted = p;
        }
    }
    /*  Compute P(barcode|observed)
        P(b|r), the probability that b was sequenced given r was observed
        P(b|r) = P(r|b) * P(b) / ( P(noise) * P(r|noise) + sigma )
        where sigma = sum of P(r|b) * P(b) over b */
    decoding_probability = adjusted / (sigma + adjusted_noise_probability);

    /* Check for decoding failure and assign to the unclassified channel if decoding failed */
    if(!(conditioned_decoding_probability > random_barcode_probability && decoding_probability > confidence_threshold)) {
        this->decoding_distance = 0;
        this->decoding_probability = 0;
        this->decoded = &this->unclassified;
    }
};

MultiplexMDDecoder::MultiplexMDDecoder(const Value& ontology) try :
    MDDecoder< Channel >(ontology) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("MultiplexMDDecoder :: " + error.message);

    } catch(exception& error) {
        throw InternalError("MultiplexMDDecoder :: " + string(error.what()));
};
void MultiplexMDDecoder::decode(const Read& input, Read& output) {
    MDDecoder< Channel >::decode(input, output);
    output.assign_RG(this->decoded->rg);
    output.update_multiplex_barcode(this->observation);
    output.update_multiplex_distance(this->decoding_distance);
};

MultiplexPAMLDecoder::MultiplexPAMLDecoder(const Value& ontology) try :
    PAMLDecoder< Channel >(ontology) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("MultiplexPAMLDecoder :: " + error.message);

    } catch(exception& error) {
        throw InternalError("MultiplexPAMLDecoder :: " + string(error.what()));
};
void MultiplexPAMLDecoder::decode(const Read& input, Read& output) {
    PAMLDecoder< Channel >::decode(input, output);
    output.assign_RG(this->decoded->rg);
    output.update_multiplex_barcode(this->observation);
    output.update_multiplex_distance(this->decoding_distance);
    output.update_multiplex_decoding_confidence(this->decoding_probability);
};

CellularMDDecoder::CellularMDDecoder(const Value& ontology) try :
    MDDecoder< Barcode >(ontology) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("CellularMDDecoder :: " + error.message);

    } catch(exception& error) {
        throw InternalError("CellularMDDecoder :: " + string(error.what()));
};
void CellularMDDecoder::decode(const Read& input, Read& output) {
    MDDecoder< Barcode >::decode(input, output);
    output.update_cellular_barcode(*this->decoded);
    output.update_raw_cellular_barcode(this->observation);
    if(this->decoded != &this->unclassified) {
        output.update_cellular_distance(this->decoding_distance);
    } else {
        output.set_cellular_distance(0);
    }
};

CellularPAMLDecoder::CellularPAMLDecoder(const Value& ontology) try :
    PAMLDecoder< Barcode >(ontology) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("CellularPAMLDecoder :: " + error.message);

    } catch(exception& error) {
        throw InternalError("CellularPAMLDecoder :: " + string(error.what()));
};
void CellularPAMLDecoder::decode(const Read& input, Read& output) {
    PAMLDecoder< Barcode >::decode(input, output);
    output.update_raw_cellular_barcode(this->observation);
    output.update_cellular_barcode(*this->decoded);
    if(this->decoded != &this->unclassified) {
        output.update_cellular_decoding_confidence(this->decoding_probability);
        output.update_cellular_distance(this->decoding_distance);
    } else {
        output.set_cellular_decoding_confidence(0);
        output.set_cellular_distance(0);
    }
};

MolecularNaiveDecoder::MolecularNaiveDecoder(const Value& ontology) try :
    Decoder(ontology),
    nucleotide_cardinality(decode_value_by_key< int32_t >("nucleotide cardinality", ontology)),
    rule(decode_value_by_key< Rule >("transform", ontology)),
    observation(decode_value_by_key< int32_t >("segment cardinality", ontology)) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("MolecularNaiveDecoder :: " + error.message);

    } catch(exception& error) {
        throw InternalError("MolecularNaiveDecoder :: " + string(error.what()));
};
void MolecularNaiveDecoder::decode(const Read& input, Read& output) {
    observation.clear();
    rule.apply(input, observation);
    output.update_molecular_barcode(observation);
};
