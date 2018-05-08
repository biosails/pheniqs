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

#ifndef PHENIQS_JSON_H
#define PHENIQS_JSON_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>

#define RAPIDJSON_NO_SIZETYPEDEFINE
namespace rapidjson { typedef ::std::size_t SizeType; }
#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/error/en.h>

#include "error.h"
#include "kstring.h"

using std::string;
using std::vector;
using std::list;
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
using rapidjson::StringRef;
using rapidjson::kNullType;
using rapidjson::kStringType;
using rapidjson::kArrayType;
using rapidjson::kObjectType;

void print_json(const Value& node, ostream& o);

/*  Recursively merge two JSON documents.
    Overlay < other > on < node > and write the result to < container >, using < document > for memory allocation.
    < container > is first initialized to a null value.
*/
void merge_json_value(const Value& node, const Value& other, Value& container, Document& document);
void merge_json_value(const Value& node, const Value& other, Document& document);

/* decoding JSON container into an object */
template < typename T > bool decode_value(T& value, const Value& container);
template < typename T > bool decode_value_by_key(const Value::Ch* key, T& value, const Value& container);
template < typename T > T decode_value_by_key(const Value::Ch* key, const Value& container) {
    T value;
    decode_value_by_key(key, value, container);
    return value;
};

/* encoding object to JSON container */
inline bool encode_key_value(const string& key, const bool& value, Value& container, Document& document) {
    container.RemoveMember(key.c_str());
    container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(value).Move(), document.GetAllocator());
    return true;
};
inline bool encode_key_value(const string& key, const uint8_t& value, Value& container, Document& document) {
    if(value < numeric_limits< uint8_t >::max()) {
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(static_cast< unsigned>(value)).Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const vector< uint8_t >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(static_cast< unsigned>(v)).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const list< uint8_t >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(static_cast< unsigned>(v)).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const int32_t& value, Value& container, Document& document) {
    if(value < numeric_limits< int32_t >::max()) {
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(static_cast< int >(value)).Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const vector< int32_t >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const uint32_t& value, Value& container, Document& document) {
    if(value < numeric_limits< uint32_t >::max()) {
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(static_cast< unsigned>(value)).Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const vector< uint32_t >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(static_cast< unsigned>(v)).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const int64_t& value, Value& container, Document& document) {
    if(value < numeric_limits< int64_t >::max()) {
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(value).Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const vector< int64_t >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const list< int64_t >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const uint64_t& value, Value& container, Document& document) {
    if(value < numeric_limits< uint64_t >::max()) {
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(value).Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const vector< uint64_t >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const list< uint64_t >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const double& value, Value& container, Document& document) {
    if(value != numeric_limits< double >::infinity()) {
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(value).Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const vector< double >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const list< double >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const string& value, Value& container, Document& document) {
    if(!value.empty()) {
        container.RemoveMember(key.c_str());
        container.AddMember(
            Value(key.c_str(), key.size(), document.GetAllocator()).Move(),
            Value(value.c_str(),value.length(), document.GetAllocator()).Move(),
            document.GetAllocator()
        );
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const vector< string >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v.c_str(), v.length(), document.GetAllocator()).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const list< string >& value, Value& container, Document& document) {
    if(!value.empty()) {
        Value array(kArrayType);
        for(auto& v : value) {
            array.PushBack(Value(v.c_str(), v.length(), document.GetAllocator()).Move(), document.GetAllocator());
        }
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
        return true;
    }
    return false;
};
inline bool encode_key_value(const string& key, const kstring_t& value, Value& container, Document& document) {
    if(value.l > 0) {
        container.RemoveMember(key.c_str());
        container.AddMember(
            Value(key.c_str(), key.size(), document.GetAllocator()).Move(),
            Value(value.s, value.l, document.GetAllocator()).Move(),
            document.GetAllocator()
        );
        return true;
    }
    return false;
};

/* transcoding from one JSON container to another */
template < typename T > bool transcode_value(const Value& from, Value& to, Document& document) {
    T buffer;
    if(decode_value< T >(buffer, from)) {
        return encode_value(buffer, to, document);
    }
    return false;
};
template < typename T > bool transcode_value_by_key(const Value::Ch* key, const Value& from, Value& to, Document& document) {
    T buffer;
    if(decode_value_by_key< T >(key, buffer, from)) {
        return encode_key_value(key, buffer, to, document);
    }
    return false;
};

#endif /* PHENIQS_JSON_H */
