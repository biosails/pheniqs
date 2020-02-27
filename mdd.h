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

#ifndef PHENIQS_MDD_H
#define PHENIQS_MDD_H

#include "include.h"
#include "decoder.h"

template < class T > class MdDecoder : public ObservingDecoder< T > {
    protected:
        const uint8_t quality_masking_threshold;
        const vector< int32_t > distance_tolerance;
        unordered_map< string, T* > element_by_sequence;

    public:
        MdDecoder(const Value& ontology);
        inline void classify(const Read& input, Read& output) override;

    private:
        inline bool match(T& barcode);
};

class MdSampleDecoder : public MdDecoder< Barcode > {
    public:
        vector< string > rg_by_barcode_index;
        MdSampleDecoder(const Value& ontology);
        inline void classify(const Read& input, Read& output) override;
};

class MdCellularDecoder : public MdDecoder< Barcode > {
    public:
        MdCellularDecoder(const Value& ontology);
        inline void classify(const Read& input, Read& output) override;
};

class MdMolecularDecoder : public MdDecoder< Barcode > {
    public:
        MdMolecularDecoder(const Value& ontology);
        inline void classify(const Read& input, Read& output) override;
};
#endif /* PHENIQS_MDD_H */
