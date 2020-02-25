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
