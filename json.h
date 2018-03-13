/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2017  Lior Galanti
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

#ifndef PHENIQS_JSON_H
#define PHENIQS_JSON_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "error.h"

#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/error/en.h>

using std::size_t;
using std::string;
using std::vector;
using std::ostream;
using std::ifstream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::numeric_limits;
using std::unordered_set;

using rapidjson::Document;
using rapidjson::Value;
using rapidjson::SizeType;
using rapidjson::StringBuffer;
using rapidjson::PrettyWriter;

void merge_json_value(const Value& node, const Value& other, Value& container, Document& document);
void merge_json_value(const Value& node, const Value& other, Document& document);

/*
    JSON decoding
*/
inline string* decode_string_by_key(const Value::Ch* key, const Value& container) {
    string* value = NULL;
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsString()) {
            value = new string();
            value->assign(element->value.GetString(), element->value.GetStringLength());
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return value;
};
inline void decode_string_by_key(const Value::Ch* key, string& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsString()) {
            value.assign(element->value.GetString(), element->value.GetStringLength());
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
};
inline void decode_uint8_by_key(const Value::Ch* key, uint8_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsUint()) {
            uint32_t v = element->value.GetUint();
            if(v < numeric_limits< uint8_t >::max()) {
                value = v;
            } else { throw ConfigurationError(string(key) + " element must be an integer smaller than " + to_string(numeric_limits< uint8_t >::max())); }
        } else { throw ConfigurationError(string(key) + " element must be an unsigned 32 bit integer"); }
    }
};
inline void decode_uint32_by_key(const Value::Ch* key, int64_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsUint()) {
            value = element->value.GetUint();
        } else { throw ConfigurationError(string(key) + " element must be an unsigned 32 bit integer"); }
    }
};
inline void decode_uint32_by_key(const Value::Ch* key, uint32_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsUint()) {
            value = element->value.GetUint();
        } else { throw ConfigurationError(string(key) + " element must be an unsigned 32 bit integer"); }
    }
};
inline void decode_uint32_by_key(const Value::Ch* key, int32_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsUint()) {
            uint32_t v = element->value.GetUint();
            if(v < numeric_limits< int32_t >::max()) {
                value = v;
            } else { throw ConfigurationError(string(key) + " element must be an integer smaller than " + to_string(numeric_limits< int32_t >::max())); }
        } else { throw ConfigurationError(string(key) + " element must be an unsigned 32 bit integer"); }
    }
};
inline void decode_int32_by_key(const Value::Ch* key, int32_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsInt()) {
            value = element->value.GetInt();
        } else { throw ConfigurationError(string(key) + " element must be a 32 bit integer"); }
    }
};
inline void decode_int32_by_key(const Value::Ch* key, int64_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsInt()) {
            value = element->value.GetInt();
        } else { throw ConfigurationError(string(key) + " element must be a 32 bit integer"); }
    }
};
inline void decode_int64_by_key(const Value::Ch* key, int64_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsInt64()) {
            value = element->value.GetInt64();
        } else { throw ConfigurationError(string(key) + " element must be a 64 bit integer"); }
    }
};
inline void decode_bool_by_key(const Value::Ch* key, bool& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsBool()) {
            value = element->value.GetBool();
        } else { throw ConfigurationError(string(key) + " element must be a boolean"); }
    }
};
inline void decode_double_by_key(const Value::Ch* key, double& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsNumber()) {
            value = element->value.GetDouble();
        } else { throw ConfigurationError(string(key) + " element must be numeric"); }
    }
};
inline void decode_double_vector_by_key(const Value::Ch* key, vector< double >& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsArray()) {
            for(const auto& v : element->value.GetArray()) {
                if(v.IsNumber()) {
                    value.push_back(v.GetDouble());
                } else { throw ConfigurationError(string(key) + " element must be numeric"); }
            }
        }
    }
};
inline void decode_size_t_by_key(const Value::Ch* key, size_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd()) {
        if(element->value.IsUint64()) {
            uint64_t v = element->value.GetUint64();
            if(v < numeric_limits< size_t >::max()) {
                value = v;
            } else { throw ConfigurationError(string(key) + " element must be an integer smaller than " + to_string(numeric_limits< size_t >::max())); }
        } else { throw ConfigurationError(string(key) + " element must be an unsigned 64 bit integer"); }
    }
};

/*
    JSON encoding
*/
inline void encode_key_value(const string& key, const double& value, Value& container, Document& document) {
    if(value != numeric_limits< double >::infinity()) {
        Value v(value);
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.AddMember(k.Move(), v.Move(), document.GetAllocator());
    }
};
inline void encode_key_value(const string& key, const uint32_t& value, Value& container, Document& document) {
    if(value < numeric_limits< int32_t >::max()) {
        Value v(value);
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.AddMember(k.Move(), v.Move(), document.GetAllocator());
    }
};
inline void encode_key_value(const string& key, const int64_t& value, Value& container, Document& document) {
    if(value < numeric_limits< int64_t >::max()) {
        Value v(value);
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.AddMember(k.Move(), v.Move(), document.GetAllocator());
    }
};
inline void encode_key_value(const string& key, const int32_t& value, Value& container, Document& document) {
    if(value < numeric_limits< int32_t >::max()) {
        Value v(value);
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.AddMember(k.Move(), v.Move(), document.GetAllocator());
    }
};
inline void encode_key_value(const string& key, const size_t& value, Value& container, Document& document) {
    if(value < numeric_limits< size_t >::max()) {
        Value v;
        v.SetUint64(value);
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.AddMember(k.Move(), v.Move(), document.GetAllocator());
    }
};
inline void encode_key_value(const string& key, const string& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value v(value.c_str(), value.length(), document.GetAllocator());
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.AddMember(k.Move(), v.Move(), document.GetAllocator());
    }
};
inline void encode_key_value(const string& key, const bool& value, Value& container, Document& document) {
    Value v(value);
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
inline void encode_key_value(const string& key, const vector< double >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value double_array;
        double_array.SetArray();
        for(auto& v : value) {
            double_array.PushBack(v, document.GetAllocator());
        }
        Value k(key.c_str(), key.size(), document.GetAllocator());
        container.AddMember(k.Move(), double_array.Move(), document.GetAllocator());
    }
};

#endif /* PHENIQS_JSON_H */