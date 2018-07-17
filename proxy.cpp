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

#include "proxy.h"

void to_string(const FormatKind& value, string& result) {
    switch(value) {
        case FormatKind::FASTQ:         result.assign("FASTQ");      break;
        case FormatKind::HTS:           result.assign("HTS");        break;
        case FormatKind::DEV_NULL:      result.assign("DEV_NULL");   break;
        default:                        result.assign("UNKNOWN");    break;
    }
};
bool from_string(const char* value, FormatKind& result) {
         if(value == NULL)              result = FormatKind::UNKNOWN;
    else if(!strcmp(value, "FASTQ"))    result = FormatKind::FASTQ;
    else if(!strcmp(value, "HTS"))      result = FormatKind::HTS;
    else if(!strcmp(value, "DEV_NULL"))      result = FormatKind::DEV_NULL;
    else                                result = FormatKind::UNKNOWN;

    return (result == FormatKind::UNKNOWN ? false : true);
};
void to_kstring(const FormatKind& value, kstring_t& result) {
    ks_clear(result);
    string string_value;
    to_string(value, string_value);
    ks_put_string(string_value.c_str(), string_value.size(), result);
};
bool from_string(const string& value, FormatKind& result) {
    return from_string(value.c_str(), result);
};
ostream& operator<<(ostream& o, const FormatKind& value) {
    string string_value;
    to_string(value, string_value);
    o << string_value;
    return o;
};
void encode_key_value(const string& key, const FormatKind& value, Value& container, Document& document) {
    string string_value;
    to_string(value, string_value);
    Value v(string_value.c_str(), string_value.length(), document.GetAllocator());
    Value k(key.c_str(), key.size(), document.GetAllocator());
    container.RemoveMember(key.c_str());
    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
};
template<> bool decode_value_by_key< FormatKind >(const Value::Ch* key, FormatKind& value, const Value& container) {
    Value::ConstMemberIterator element = container.FindMember(key);
    if(element != container.MemberEnd() && !element->value.IsNull()) {
        if(element->value.IsString()) {
            return from_string(element->value.GetString(), value);
        } else { throw ConfigurationError(string(key) + " element must be a string"); }
    }
    return false;
};

bool encode_value(const string& key, const FeedProxy& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        Value element(kObjectType);
        encode_key_value("index", value.index, element, document);
        encode_key_value("url", value.url, element, document);
        encode_key_value("direction", value.direction, element, document);
        encode_key_value("platform", value.platform, element, document);
        encode_key_value("capacity", value.capacity, element, document);
        encode_key_value("resolution", value.resolution, element, document);
        encode_key_value("phred offset", value.phred_offset, element, document);
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), element.Move(), document.GetAllocator());
        return true;
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
ostream& operator<<(ostream& o, const FeedProxy& proxy) {
    o << "direction : " << proxy.direction << endl;
    o << "index : " << proxy.index << endl;
    o << "url : " << proxy.url << endl;
    o << "platform : " << proxy.platform << endl;
    o << "capacity : " << proxy.capacity << endl;
    o << "resolution : " << proxy.resolution << endl;
    o << "phred_offset : " << to_string(proxy.phred_offset) << endl;
    proxy.url.describe(o);
    return o;
};
template<> FeedProxy decode_value< FeedProxy >(const Value& container) {
    if(container.IsObject()) {
        FeedProxy proxy(container);
        return proxy;
    } else { throw ConfigurationError("feed proxy element must be a dictionary"); }
};
template<> list< FeedProxy > decode_value_by_key< list< FeedProxy > >(const Value::Ch* key, const Value& container) {
    if(container.IsObject()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            if(!reference->value.IsNull()) {
                if(reference->value.IsArray()) {
                    list< FeedProxy > value;
                    if(!reference->value.Empty()) {
                        for(auto& element : reference->value.GetArray()) {
                            value.emplace_back(element);
                        }
                    }
                    return value;
                } else { throw ConfigurationError(string(key) + " is not an array"); }
            } else { throw ConfigurationError(string(key) + " is null"); }
        } else { throw ConfigurationError(string(key) + " not found"); }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
};
