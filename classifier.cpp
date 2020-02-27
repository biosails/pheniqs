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

#include "classifier.h"

/* ClassifierType */

string to_string(const ClassifierType& value) {
    string result;
    switch(value) {
    case ClassifierType::SAMPLE:      result.assign("sample");    break;
    case ClassifierType::CELLULAR:    result.assign("cellular");  break;
    case ClassifierType::MOLECULAR:   result.assign("molecular"); break;
        default:                        result.assign("unknown");   break;
    }
    return result;
};
bool from_string(const char* value, ClassifierType& result) {
         if(value == NULL)                  result = ClassifierType::UNKNOWN;
    else if(!strcmp(value, "sample"))       result = ClassifierType::SAMPLE;
    else if(!strcmp(value, "cellular"))     result = ClassifierType::CELLULAR;
    else if(!strcmp(value, "molecular"))    result = ClassifierType::MOLECULAR;
    else                                    result = ClassifierType::UNKNOWN;

    return (result == ClassifierType::UNKNOWN ? false : true);
};
bool from_string(const string& value, ClassifierType& result) {
    return from_string(value.c_str(), result);
};
void to_kstring(const ClassifierType& value, kstring_t& result) {
    ks_clear(result);
    string string_value(to_string(value));
    ks_put_string(string_value.c_str(), string_value.size(), result);
};
ostream& operator<<(ostream& o, const ClassifierType& value) {
    o << to_string(value);
    return o;
};
void encode_key_value(const string& key, const ClassifierType& value, Value& container, Document& document) {
    string string_value(to_string(value));
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< ClassifierType >(const Value::Ch* key, ClassifierType& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
template <> ClassifierType decode_value_by_key(const Value::Ch* key, const Value& container) {
    ClassifierType value(ClassifierType::UNKNOWN);
    decode_value_by_key(key, value, container);
    return value;
};

vector< string > decode_tag_ID_by_index(const Value& ontology) {
    vector< string > value;
    if(ontology.IsObject()) {
        Value::ConstMemberIterator undetermined_reference = ontology.FindMember("undetermined");
        if(undetermined_reference != ontology.MemberEnd()) {
            Value::ConstMemberIterator codec_reference = ontology.FindMember("codec");
            if(codec_reference != ontology.MemberEnd()) {
                value.reserve(codec_reference->value.MemberCount() + 1);
                value.emplace_back(decode_value_by_key< string >("ID", undetermined_reference->value));
                for(auto& record : codec_reference->value.GetObject()) {
                    value.emplace_back(decode_value_by_key< string >("ID", record.value));
                }
            } else {
                value.reserve(1);
                value.emplace_back(decode_value_by_key< string >("ID", undetermined_reference->value));
            }
        } else { throw ConfigurationError("classifier must declare an undetermined element"); }
    }
    return value;
};
