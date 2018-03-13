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

#include "json.h"

void merge_json_value(const Value& node, const Value& other, Value& container, Document& document) {
    container.SetNull();
    if(!other.IsNull() && !node.IsNull()) {
        if(other.IsObject() && node.IsObject()) {
            container.SetObject();
            for(auto& record : node.GetObject()) {
                Value::ConstMemberIterator element = other.FindMember(record.name);
                if(element != other.MemberEnd()) {
                    Value next;
                    merge_json_value(record.value, element->value, next, document);
                    Value k(record.name, document.GetAllocator());
                    container.AddMember(k.Move(), next.Move(), document.GetAllocator());
                } else {
                    Value v;
                    v.CopyFrom(record.value, document.GetAllocator());
                    Value k(record.name, document.GetAllocator());
                    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
                }
            }
            for(auto& record : other.GetObject()) {
                Value::ConstMemberIterator element = container.FindMember(record.name);
                if(element == container.MemberEnd()) {
                    Value v;
                    v.CopyFrom(record.value, document.GetAllocator());
                    Value k(record.name, document.GetAllocator());
                    container.AddMember(k.Move(), v.Move(), document.GetAllocator());
                }
            }
        } else if(other.IsArray() && node.IsArray()) {
            container.CopyFrom(node, document.GetAllocator());
            for(const auto& e : other.GetArray()) {
                Value v;
                v.CopyFrom(e, document.GetAllocator());
                container.PushBack(v.Move(), document.GetAllocator());
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
