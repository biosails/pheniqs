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

/* Molecular */
NaiveMolecularDecoder::NaiveMolecularDecoder(const Value& ontology) try :
    ObservingDecoder< Barcode >(ontology) {

    } catch(Error& error) {
        error.push("NaiveMolecularDecoder");
        throw;
};
void NaiveMolecularDecoder::decode(const Read& input, Read& output) {
    this->observation.clear();
    this->rule.apply(input, this->observation);
    output.update_molecular_barcode(observation);
    ObservingDecoder< Barcode >::decode(input, output);
};
