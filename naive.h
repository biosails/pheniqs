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

#ifndef PHENIQS_NAIVE_H
#define PHENIQS_NAIVE_H

#include "include.h"
#include "decoder.h"

class NaiveMolecularDecoder : public Decoder< Barcode > {
    public:
        NaiveMolecularDecoder(const Value& ontology) try :
            Decoder< Barcode >(ontology) {

            } catch(Error& error) {
                error.push("NaiveMolecularDecoder");
                throw;
        };
        NaiveMolecularDecoder(const NaiveMolecularDecoder& other) :
            Decoder< Barcode >(other) {
        };
        inline void classify(const Read& input, Read& output) override {
            this->observation.clear();
            this->rule.apply(input, this->observation);
            output.update_molecular_barcode(observation);
            Decoder< Barcode >::classify(input, output);
        };
};

#endif /* PHENIQS_NAIVE_H */
