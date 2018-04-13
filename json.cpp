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

#include "json.h"


void print_json(const Value& node, ostream& o) {
    StringBuffer buffer;
    PrettyWriter< StringBuffer > writer(buffer);
    node.Accept(writer);
    o << buffer.GetString() << endl;
};

/*  Recursively merge two JSON documents.
    This will overlay other on node and write the result to container, using document for memory allocation.
*/
void merge_json_value(const Value& node, const Value& other, Value& container, Document& document) {
    container.SetNull();
    if(!other.IsNull() && !node.IsNull()) {
        if(other.IsObject() && node.IsObject()) {
            container.SetObject();
            for(auto& record : node.GetObject()) {
                Value::ConstMemberIterator element = other.FindMember(record.name);
                Value next;
                if(element != other.MemberEnd()) {
                    merge_json_value(record.value, element->value, next, document);
                } else {
                    next.CopyFrom(record.value, document.GetAllocator());
                }
                container.AddMember(Value(record.name, document.GetAllocator()).Move(), next.Move(), document.GetAllocator());
            }
            for(auto& record : other.GetObject()) {
                Value::ConstMemberIterator element = container.FindMember(record.name);
                if(element == container.MemberEnd()) {
                    Value next;
                    next.CopyFrom(record.value, document.GetAllocator());
                    container.AddMember(Value(record.name, document.GetAllocator()).Move(), next.Move(), document.GetAllocator());
                }
            }
        } else if(other.IsArray() && node.IsArray()) {
            container.CopyFrom(node, document.GetAllocator());
            for(const auto& e : other.GetArray()) {
                Value next;
                next.CopyFrom(e, document.GetAllocator());
                container.PushBack(next.Move(), document.GetAllocator());
            }
        }
    }
    if(!other.IsNull() && container.IsNull()) {
        container.CopyFrom(other, document.GetAllocator());
    }
};
void merge_json_value(const Value& node, const Value& other, Document& document) {
    merge_json_value(node, other, document, document);
};

template<> bool decode_value_by_key< bool >(const Value::Ch* key, bool& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsBool()) {
            value = element->value.GetBool();
            return true;
        } else { throw ConfigurationError(string(key) + " element must be a boolean"); }
    }
    return false;
};
template<> bool decode_value_by_key< uint8_t >(const Value::Ch* key, uint8_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        unsigned v;
        if(element->value.IsUint() && (v = element->value.GetUint()) < numeric_limits< uint8_t >::max()) {
            value = static_cast< uint8_t >(v);
            return true;
        } else { throw ConfigurationError(string(key) + " element must be an unsigned integer smaller than " + to_string(numeric_limits< uint8_t >::max())); }
    }
    return false;
};
template<> bool decode_value_by_key< vector< uint8_t > >(const Value::Ch* key, vector< uint8_t >& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsArray()) {
            if(!element->value.Empty()) {
                value.reserve(value.size() + element->value.Size());
                unsigned v;
                for(const auto& e : element->value.GetArray()) {
                    if(e.IsUint() && (v = e.GetUint()) < numeric_limits< uint8_t >::max()) {
                        value.emplace_back(static_cast < uint8_t >(v));
                    } else { throw ConfigurationError(string(key) + " element member must be an unsigned integer smaller than " + to_string(numeric_limits< uint8_t >::max())); }
                }
                return true;
            }
        } else { throw ConfigurationError(string(key) + " element must be an array"); }
    }
    return false;
};
template<> bool decode_value_by_key< int32_t >(const Value::Ch* key, int32_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        int64_t v;
        if(element->value.IsInt64() && (v = element->value.GetInt64()) < numeric_limits< int32_t >::max() && v > numeric_limits< int32_t >::min()) {
            value = static_cast< int32_t >(v);
            return true;
        } else { throw ConfigurationError(string(key) + " element must be an integer between " + to_string(numeric_limits< int32_t >::min()) + " and " + to_string(numeric_limits< int32_t >::max())); }
    }
    return false;
};
template<> bool decode_value_by_key< uint32_t >(const Value::Ch* key, uint32_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        uint64_t v;
        if(element->value.IsUint64() && (v = element->value.GetUint64()) < numeric_limits< uint32_t >::max()) {
            value = static_cast< uint32_t >(v);
            return true;
        } else { throw ConfigurationError(string(key) + " element must be an unsigned integer smaller than " + to_string(numeric_limits< uint32_t >::max())); }
    }
    return false;
};
template<> bool decode_value_by_key< int64_t >(const Value::Ch* key, int64_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        int64_t v;
        if(element->value.IsInt64() && (v = element->value.GetInt64()) < numeric_limits< int64_t >::max()) {
            value = v;
            return true;
        } else { throw ConfigurationError(string(key) + " element must be an integer between " + to_string(numeric_limits< int64_t >::min()) + " and " + to_string(numeric_limits< int64_t >::max())); }
    }
    return false;
};
template<> bool decode_value_by_key< uint64_t >(const Value::Ch* key, uint64_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        uint64_t v;
        if(element->value.IsUint64() && (v = element->value.GetUint64()) < numeric_limits< uint64_t >::max()) {
            value = v;
            return true;
        } else { throw ConfigurationError(string(key) + " element must be an unsigned integer smaller than " + to_string(numeric_limits< uint64_t >::max())); }
    }
    return false;
};
template<> bool decode_value_by_key< vector< uint64_t > >(const Value::Ch* key, vector< uint64_t >& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsArray()){
            if(!element->value.Empty()) {
                value.reserve(value.size() + element->value.Size());
                uint64_t v;
                for(const auto& e : element->value.GetArray()) {
                    if(e.IsUint64() && (v = e.GetUint64()) < numeric_limits< uint64_t >::max()) {
                        value.emplace_back(v);
                    } else { throw ConfigurationError(string(key) + " element member must be an unsigned integer smaller than " + to_string(numeric_limits< uint64_t >::max())); }
                }
                return true;
            }
        } else { throw ConfigurationError(string(key) + " element must be an array"); }
    }
    return false;
};
template<> bool decode_value_by_key< double >(const Value::Ch* key, double& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsNumber()) {
            value = element->value.GetDouble();
            return true;
        } else { throw ConfigurationError(string(key) + " element must be numeric"); }
    }
    return false;
};
template<> bool decode_value_by_key< vector< double > >(const Value::Ch* key, vector< double >& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsArray()) {
            if(!element->value.Empty()) {
                value.reserve(value.size() + element->value.Size());
                for(const auto& e : element->value.GetArray()) {
                    if(e.IsNumber()) {
                        value.emplace_back(e.GetDouble());
                    } else { throw ConfigurationError(string(key) + " element member must be numeric"); }
                }
                return true;
            }
        } else { throw ConfigurationError(string(key) + " element must be an array"); }
    }
    return false;
};
template<> bool decode_value_by_key< string >(const Value::Ch* key, string& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            value.assign(element->value.GetString(), element->value.GetStringLength());
            return true;
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};
template<> bool decode_value_by_key< kstring_t >(const Value::Ch* key, kstring_t& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            if(element->value.GetStringLength() < numeric_limits< size_t >::max()) {
                ks_put_string(element->value.GetString(), static_cast< size_t >(element->value.GetStringLength()), value);
                return true;
            } else { throw ConfigurationError(string(key) + " element must be a string shorter than " + to_string(numeric_limits< int >::max())); }
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};

template <> bool decode_value_by_key(const Value::Ch* key, const Value& container) {
    bool value(false);
    decode_value_by_key(key, value, container);
    return value;
};
template <> uint64_t decode_value_by_key(const Value::Ch* key, const Value& container) {
    uint64_t value(numeric_limits< uint64_t >::max());
    decode_value_by_key(key, value, container);
    return value;
};
template <> double decode_value_by_key(const Value::Ch* key, const Value& container) {
    double value(numeric_limits< double >::infinity());
    decode_value_by_key(key, value, container);
    return value;
};
