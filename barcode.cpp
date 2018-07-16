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

Barcode::Barcode(const Value& ontology) :
    SequenceArray< Sequence >(decode_value_by_key< int32_t >("segment cardinality", ontology)),
    index(decode_value_by_key< int32_t >("index", ontology)),
    concentration(decode_value_by_key< double >("concentration", ontology)) {
    Value::ConstMemberIterator reference = ontology.FindMember("barcode");
    if(reference != ontology.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsArray()) {
                if(reference->value.Size() == segment_cardinality()) {
                    size_t segment_index(0);
                    for(auto& segment : reference->value.GetArray()) {
                        if(!segment.IsNull()) {
                            if(segment.IsString()) {
                                // if(is_iupac_strict_sequence(segment.GetString())) {
                                    segment_array[segment_index].fill(segment.GetString(), static_cast< int32_t >(segment.GetStringLength()));
                                    ++segment_index;
                                // } else { throw ConfigurationError("barcode segment " + to_string(segment_index) + " " + string(segment.GetString(), segment.GetStringLength()) + " is not strict IUPAC encoded"); }
                            } else { throw ConfigurationError("barcode segment " + to_string(segment_index) + " must be a string"); }
                        } else { throw ConfigurationError("barcode segment " + to_string(segment_index) + " can not be null"); }
                    }
                } else { throw ConfigurationError("barcode must have exactly " + to_string(segment_cardinality()) + " segments"); }
            } else { throw ConfigurationError("barcode segment element must be an array"); }
        } else { throw ConfigurationError("barcode segment can not be null"); }
    } else { throw ConfigurationError("barcode is missing a barcode segment element"); }
};

template<> vector< Barcode > decode_value_by_key(const Value::Ch* key, const Value& container) {
    vector< Barcode > value;
    Value::ConstMemberIterator reference = container.FindMember(key);
    if(reference != container.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                value.reserve(reference->value.MemberCount());
                for(auto& record : reference->value.GetObject()) {
                    value.emplace_back(record.value);
                }
                value.shrink_to_fit();
            } else { throw ConfigurationError(string(key) + " element must be a dictionary"); }
        } else { throw ConfigurationError(string(key) + " element is null"); }
    } else { throw ConfigurationError(string(key) + " not found"); }
    return value;
};

bool encode_key_value(const string& key, const Barcode& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value element(kObjectType);
        encode_key_value("index", value.index, element, document);
        encode_key_value("segment cardinality", static_cast< int32_t >(value.segment_cardinality()), element, document);
        encode_key_value("concentration", value.concentration, element, document);
        Value array(kArrayType);
        for(auto& segment : value) {
            encode_value(segment, array, document);
        }
        element.AddMember(Value("barcode", document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
ostream& operator<<(ostream& o, const Barcode& barcode) {
    o << barcode.iupac_ambiguity() << endl;
    return o;
};
