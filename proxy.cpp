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
