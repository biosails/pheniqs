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

#ifndef PHENIQS_CHANNEL_H
#define PHENIQS_CHANNEL_H

#include "include.h"

#include "feed.h"

class Channel : public Barcode {
    void operator=(Channel const &) = delete;

    public:
        const HeadRGAtom rg;
        const bool concrete;
        const bool disable_quality_control;
        const list< URL > output_feed_url_by_segment;
        vector< Feed* > output_feed_by_order;
        vector< Feed* > output_feed_by_segment;
        Channel(const Value& ontology) :
            Barcode(ontology),
            rg(ontology),
            concrete(!decode_value_by_key< bool >("output disabled", ontology)),
            disable_quality_control(decode_value_by_key< bool >("disable quality control", ontology)),
            output_feed_url_by_segment(decode_value_by_key< list< URL > >("output", ontology)) {
        };
        Channel(const Channel& other) :
            Barcode(other),
            rg(other.rg),
            concrete(other.concrete),
            disable_quality_control(other.disable_quality_control),
            output_feed_url_by_segment(other.output_feed_url_by_segment),
            output_feed_by_order(other.output_feed_by_order),
            output_feed_by_segment(other.output_feed_by_segment) {
        };
        inline void push(const Read& read) {
            if(concrete) {
                // acquire a push lock for all feeds in a fixed order
                vector< unique_lock< mutex > > feed_locks;
                feed_locks.reserve(output_feed_by_order.size());
                for(const auto feed : output_feed_by_order) {
                    feed_locks.push_back(feed->acquire_push_lock());
                }

                // push the segments to the output feeds
                for(size_t i = 0; i < output_feed_by_segment.size(); ++i) {
                    output_feed_by_segment[i]->push(read[i]);
                }

                // release the locks on the feeds in reverse order
                for(auto feed_lock = feed_locks.rbegin(); feed_lock != feed_locks.rend(); ++feed_lock) {
                    feed_lock->unlock();
                }
            }
        };
};

template<> vector< Channel > decode_value_by_key(const Value::Ch* key, const Value& container);

#endif /* PHENIQS_CHANNEL_H */
