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

#include "barcode.h"

Barcode::Barcode(const Value& ontology) try :
    SequenceArray< Sequence >(decode_value_by_key< int32_t >("segment cardinality", ontology)),
    BarcodeAccumulator(),
    index(decode_value_by_key< int32_t >("index", ontology)),
    concentration(decode_value_by_key< double >("concentration", ontology)) {

    Value::ConstMemberIterator reference = ontology.FindMember("barcode");
    if(reference != ontology.MemberEnd()) {
        if(reference->value.Size() == segment_cardinality()) {
            size_t segment_index(0);
            for(auto& segment : reference->value.GetArray()) {
                if(segment.IsString()) {
                    segment_array[segment_index].fill(segment.GetString(), static_cast< int32_t >(segment.GetStringLength()));
                    ++segment_index;
                } else { throw ConfigurationError("barcode segment " + to_string(segment_index) + " must be a string"); }
            }
        } else { throw ConfigurationError("barcode must have exactly " + to_string(segment_cardinality()) + " segments"); }
    }

    } catch(Error& error) {
        error.push("Barcode");
        throw;
};
Barcode::Barcode(const Barcode& other) :
    SequenceArray< Sequence >(other),
    BarcodeAccumulator(other),
    index(other.index),
    concentration(other.concentration) {
};
void Barcode::encode(Value& container, Document& document) const {
    BarcodeAccumulator::encode(container, document);
    if(container.IsObject()) {
        encode_key_value("index", index, container, document);
        if(is_classified()) {
            encode_key_value("concentration", concentration, container, document);
            Value array(kArrayType);
            for(auto& segment : *this) {
                encode_value(segment, array, document);
            }
            container.AddMember(Value("barcode", document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        }
    } else { throw ConfigurationError("element must be a dictionary"); }
};
Barcode& Barcode::operator+=(const Barcode& rhs) {
    BarcodeAccumulator::operator+=(rhs);
    return *this;
};

template<> vector< Barcode > decode_value_by_key(const Value::Ch* key, const Value& container) {
    vector< Barcode > value;
    Value::ConstMemberIterator reference = container.FindMember(key);
    if(reference != container.MemberEnd()) {
        value.reserve(reference->value.MemberCount());
        for(auto& record : reference->value.GetObject()) {
            value.emplace_back(record.value);
        }
        value.shrink_to_fit();
    }
    return value;
};
ostream& operator<<(ostream& o, const Barcode& barcode) {
    o << barcode.iupac_ambiguity() << endl;
    return o;
};
