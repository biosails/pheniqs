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

#include "channel.h"

template<> vector< Channel > decode_value_by_key(const Value::Ch* key, const Value& container) {
    vector< Channel > value;
    Value::ConstMemberIterator reference = container.FindMember(key);
    if(reference != container.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                value.reserve(reference->value.MemberCount());
                for(auto& record : reference->value.GetObject()) {
                    value.emplace_back(record.value);
                }
            } else { throw ConfigurationError(string(key) + " element must be a dictionary"); }
        } else { throw ConfigurationError(string(key) + " element is null"); }
    } else { throw ConfigurationError(string(key) + " not found"); }
    return value;
};
