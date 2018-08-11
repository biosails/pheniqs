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

Channel::Channel(const Value& ontology) try :
    Barcode(ontology),
    rg(ontology),
    include_filtered(decode_value_by_key< bool >("include filtered", ontology)),
    disable_quality_control(decode_value_by_key< bool >("disable quality control", ontology)),
    output_feed_url_by_segment(decode_value_by_key< list< URL > >("output", ontology)) {

    } catch(Error& error) {
        error.push("Channel");
        throw;
};
Channel::Channel(const Channel& other) :
    Barcode(other),
    rg(other.rg),
    include_filtered(other.include_filtered),
    disable_quality_control(other.disable_quality_control),
    output_feed_url_by_segment(other.output_feed_url_by_segment),
    output_feed_lock_order(other.output_feed_lock_order),
    output_feed_by_segment(other.output_feed_by_segment) {
};
void Channel::populate(unordered_map< URL, Feed* >& feed_by_url) {
    map< int32_t, Feed* > feed_by_index;

    /* populate the output feed by segment array */
    output_feed_by_segment.reserve(output_feed_url_by_segment.size());
    for(const auto& url : output_feed_url_by_segment) {
        Feed* feed(feed_by_url[url]);
        output_feed_by_segment.emplace_back(feed);
        if(feed_by_index.count(feed->index) == 0) {
            feed_by_index.emplace(make_pair(feed->index, feed));
        }
    }
    output_feed_by_segment.shrink_to_fit();

    /* populate the output feed lock order array */
    output_feed_lock_order.reserve(feed_by_index.size());
    for(auto& record : feed_by_index) {
        /* /dev/null is not really being written to so we don't need to lock it */
        if(!record.second->is_dev_null()) {
            output_feed_lock_order.push_back(record.second);
        }
    }
    output_feed_lock_order.shrink_to_fit();
};

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
