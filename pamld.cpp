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

#include "pamld.h"

template < class T > PamlDecoder< T >::PamlDecoder(const Value& ontology) try :
    ObservingDecoder< T >(ontology),
    noise(decode_value_by_key< double >("noise", ontology)),
    confidence_threshold(decode_value_by_key< double >("confidence threshold", ontology)),
    random_barcode_probability(1.0 / double(pow(4, (this->nucleotide_cardinality)))),
    adjusted_noise_probability(noise * random_barcode_probability),
    conditional_decoding_probability(0),
    decoding_confidence(0) {

    } catch(Error& error) {
        error.push("PamlDecoder");
        throw;
};
template < class T > void PamlDecoder< T >::classify(const Read& input, Read& output) {
    this->observation.clear();
    this->rule.apply(input, this->observation);

    /*  Compute the posterior probability P( observed | expected ) for each barcode.
        Keep track of the channel that yield the maximal prior adjusted probability.
        If r is the observed sequence and b is the barcode sequence
        P(r|b) is the probability that r was observed given that b was sequenced.
        Accumulate probabilities P(b) * P(r|b), in sigma_p using the Kahan summation algorithm
        to minimize floating point drift. see https://en.wikipedia.org/wiki/Kahan_summation_algorithm.
        p is the prior adjusted conditional probability
        d is the decoding Hamming distance
        sigma_p accumulates the prior adjusted conditional probabilities to compute the decoding confidence
    */
    double p(0);
    double y(0);
    double t(0);
    int32_t d(0);
    double sigma_p(0);
    double compensation(0);
    double conditional_probability(0);
    double adjusted_conditional_decoding_probability(0);

    for(auto& barcode : this->tag_by_index) {
        /*  The conditional probability, P(r|b), is the probability of the observation r
            given b was expected.
            P(b), barcode.concentration, is the prior probability of observing b
            p is the prior adjusted conditional probability, P(b) * P(r|b),
            sigma_p is the sum of p over b  */

        barcode.compensated_decoding_probability(this->observation, conditional_probability, d);
        p = conditional_probability * barcode.concentration;
        y = p - compensation;
        t = sigma_p + y;
        compensation = (t - sigma_p) - y;
        sigma_p = t;
        if(p > adjusted_conditional_decoding_probability) {
            this->decoded = &barcode;
            this->decoding_hamming_distance = d;
            adjusted_conditional_decoding_probability = p;
            conditional_decoding_probability = conditional_probability;
        }
    }

    /* add the prior adjusted noise probability to sigma_p */
    y = adjusted_noise_probability - compensation;
    t = sigma_p + y;
    compensation = (t - sigma_p) - y;
    sigma_p = t;

    /*  P(b|r), decoding confidence, is the posterior probability.
        The probability that barcode b was sequenced given the observation r.
        adjusted_conditional_decoding_probability is the highest prior adjusted conditional probability,
        P(r|b) of all possible b */
    decoding_confidence = adjusted_conditional_decoding_probability / sigma_p;

    /*  This is a noise filter, when the conditional probability is lower than the probability of
        a random abservation, the entropy is too high for the information to be meaningful */
    if(conditional_decoding_probability > random_barcode_probability) {

        /*  if the posterior probability is higher than the confidence_threshold */
        if(decoding_confidence > confidence_threshold) {
            this->decoded->accumulated_confidence += decoding_confidence;
            if(!output.qcfail()) {
                this->decoded->accumulated_pf_confidence += decoding_confidence;
            }

        } else {
            ++this->decoded->low_confidence_count;
            output.set_qcfail(true);
        }

    } else {
        ++this->decoded->low_conditional_confidence_count;
        output.set_qcfail(true);
        this->decoded = &this->unclassified;
        this->decoding_hamming_distance = 0;
        decoding_confidence = 0;
    }
    ObservingDecoder< T >::classify(input, output);
};

PamlSampleDecoder::PamlSampleDecoder(const Value& ontology) try :
    PamlDecoder< Barcode >(ontology),
    rg_by_barcode_index(decode_tag_ID_by_index(ontology)) {

    } catch(Error& error) {
        error.push("PamlSampleDecoder");
        throw;
};
void PamlSampleDecoder::classify(const Read& input, Read& output) {
    PamlDecoder< Barcode >::classify(input, output);
    output.update_sample_barcode(this->observation);
    output.update_sample_distance(this->decoding_hamming_distance);
    output.update_sample_decoding_confidence(this->decoding_confidence);
    output.set_RG(this->rg_by_barcode_index[this->decoded->index]);
};

PamlCellularDecoder::PamlCellularDecoder(const Value& ontology) try :
    PamlDecoder< Barcode >(ontology) {

    } catch(Error& error) {
        error.push("PamlCellularDecoder");
        throw;
};
void PamlCellularDecoder::classify(const Read& input, Read& output) {
    PamlDecoder< Barcode >::classify(input, output);
    output.update_raw_cellular_barcode(this->observation);
    output.update_cellular_barcode(*this->decoded);
    if(this->decoded->is_classified()) {
        output.update_cellular_decoding_confidence(this->decoding_confidence);
        output.update_cellular_distance(this->decoding_hamming_distance);
    } else {
        output.set_cellular_decoding_confidence(0);
        output.set_cellular_distance(0);
    }
};
