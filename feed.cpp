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

#include "feed.h"

Value encode_value(const Feed& value, Document& document) {
    Value element(kObjectType);
    encode_key_value("index", value.index, element, document);
    encode_key_value("url", value.url, element, document);
    encode_key_value("direction", value.direction, element, document);
    encode_key_value("platform", value.platform, element, document);
    encode_key_value("capacity", value.capacity(), element, document);
    encode_key_value("resolution", value.resolution(), element, document);
    encode_key_value("phred offset", value.phred_offset, element, document);
    return element;
};
bool encode_key_value(const string& key, const list< Feed* >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& feed : value) {
                Value element(encode_value(*(feed), document));
                array.PushBack(element.Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
bool encode_key_value(const string& key, const vector< Feed* >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& feed : value) {
                Value element(encode_value(*(feed), document));
                array.PushBack(element.Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
template< typename T > ostream& operator<<(ostream& o, const CyclicBuffer< T >& buffer) {
    o << "Next: " << buffer._next << endl;
    o << "Vacant: " << buffer._vacant << endl;
    o << "Capacity: " << buffer._capacity << endl;
    o << "Resolution: " << buffer._resolution << endl;
    o << "Size: " << buffer.size() << endl;
    return o;
};
