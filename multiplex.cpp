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

#include "multiplex.h"

/*  Channel */
Channel::Channel(const Value& ontology) try :
    index(decode_value_by_key< int32_t >("index", ontology)),
    // rg(ontology),
    filter_outgoing_qc_fail(decode_value_by_key< bool >("filter outgoing qc fail", ontology)),
    enable_quality_control(decode_value_by_key< bool >("enable quality control", ontology)),
    output_feed_url_by_segment(decode_value_by_key< list< URL > >("output", ontology)),
    read_accumulator(decode_value_by_key< int32_t >("segment cardinality", ontology)) {

    } catch(Error& error) {
        error.push("Channel");
        throw;
};
Channel::Channel(const Channel& other) :
    index(other.index),
    // rg(other.rg),
    filter_outgoing_qc_fail(other.filter_outgoing_qc_fail),
    enable_quality_control(other.enable_quality_control),
    output_feed_url_by_segment(other.output_feed_url_by_segment),
    output_feed_lock_order(other.output_feed_lock_order),
    output_feed_by_segment(other.output_feed_by_segment),
    read_accumulator(other.read_accumulator) {
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
void Channel::finalize() {
    read_accumulator.finalize();
};
void Channel::encode(Value& container, Document& document) const {
    if(container.IsObject()) {
        // encode_value(rg, container, document);
        Value quality_control_by_segment(kArrayType);
        /*
        size_t index(0);
        for(const auto& url : output_feed_url_by_segment) {
            Value segment_report(kObjectType);
            encode_value(read_accumulator.segment_accumulator_by_index[index], segment_report, document);
            encode_key_value("url", url, segment_report, document);
            quality_control_by_segment.PushBack(segment_report.Move(), document.GetAllocator());
            ++index;
        }
        */
        encode_value(read_accumulator, quality_control_by_segment, document);
        container.AddMember("quality control by segment", quality_control_by_segment.Move(), document.GetAllocator());
    } else { throw ConfigurationError("element must be a dictionary"); }
};
Channel& Channel::operator+=(const Channel& rhs) {
    read_accumulator += rhs.read_accumulator;
    return *this;
};
template<> vector< Channel > decode_value_by_key(const Value::Ch* key, const Value& container) {
    vector< Channel > value;
    Value::ConstMemberIterator decoder_reference = container.FindMember(key);
    if(decoder_reference != container.MemberEnd()) {
        Value::ConstMemberIterator undetermined_reference = decoder_reference->value.FindMember("undetermined");
        if(undetermined_reference != decoder_reference->value.MemberEnd()) {
            Value::ConstMemberIterator codec_reference = decoder_reference->value.FindMember("codec");
            if(codec_reference != decoder_reference->value.MemberEnd()) {
                value.reserve(codec_reference->value.MemberCount() + 1);
                value.emplace_back(undetermined_reference->value);
                for(auto& record : codec_reference->value.GetObject()) {
                    value.emplace_back(record.value);
                }
            } else {
                value.reserve(1);
                value.emplace_back(undetermined_reference->value);
            }
        } else { throw ConfigurationError("decoder must declare an undetermined element"); }
    }
    return value;
};

/*  Multiplexer */
Multiplexer::Multiplexer(const Value& ontology) try :
    filter_outgoing_qc_fail(decode_value_by_key< bool >("filter outgoing qc fail", ontology)),
    enable_quality_control(decode_value_by_key< bool >("enable quality control", ontology)),
    channel_by_index(decode_value_by_key< vector< Channel > >("multiplex", ontology)) {

    } catch(Error& error) {
        error.push("Multiplexer");
        throw;
};
Multiplexer::Multiplexer(const Multiplexer& other) :
    filter_outgoing_qc_fail(other.filter_outgoing_qc_fail),
    enable_quality_control(other.enable_quality_control),
    channel_by_index(other.channel_by_index) {
};
Multiplexer& Multiplexer::operator+=(const Multiplexer& rhs) {
    if(enable_quality_control) {
        for(size_t index(0); index < channel_by_index.size(); ++index) {
            channel_by_index[index] += rhs.channel_by_index[index];
        }
    }
    return *this;
};
void Multiplexer::finalize() {
    if(enable_quality_control) {
        for(auto& channel : channel_by_index) {
            channel.finalize();
        }
    }
};
void Multiplexer::encode(Value& container, Document& document) const {
    if(enable_quality_control) {
        if(container.IsObject()) {
            Value channel_array(kArrayType);
            for(auto& channel : channel_by_index) {
                Value channel_report(kObjectType);
                channel.encode(channel_report, document);
                channel_array.PushBack(channel_report.Move(), document.GetAllocator());
            }
            container.AddMember("multiplex", channel_array.Move(), document.GetAllocator());
        } else { throw ConfigurationError("element must be a dictionary"); }
    }
};
