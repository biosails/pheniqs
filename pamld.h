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

#ifndef PHENIQS_PAMLD_H
#define PHENIQS_PAMLD_H

#include "include.h"
#include "decoder.h"

template < class T > class PamlDecoder : public ObservingDecoder< T > {
    protected:
        const double noise;
        const double confidence_threshold;
        const double random_barcode_probability;
        const double adjusted_noise_probability;
        double conditional_decoding_probability;
        double decoding_confidence;

    public:
        PamlDecoder(const Value& ontology);
        inline void classify(const Read& input, Read& output) override;
        inline void finalize() override {
            for(auto& element : this->tag_by_index) {
                this->accumulated_classified_confidence += element.accumulated_confidence;
                this->accumulated_pf_classified_confidence += element.accumulated_pf_confidence;
                this->low_conditional_confidence_count += element.low_conditional_confidence_count;
                this->low_confidence_count += element.low_confidence_count;
            }
            ObservingDecoder< T >::finalize();
        };
};

class PamlMultiplexDecoder : public PamlDecoder< Channel > {
    public:
        PamlMultiplexDecoder(const Value& ontology);
        inline void classify(const Read& input, Read& output) override;
};

class PamlCellularDecoder : public PamlDecoder< Barcode > {
    public:
        PamlCellularDecoder(const Value& ontology);
        inline void classify(const Read& input, Read& output) override;
};

#endif /* PHENIQS_PAMLD_H */
