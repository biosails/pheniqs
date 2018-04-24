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

#ifndef PHENIQS_MODEL_H
#define PHENIQS_MODEL_H

#include <set>
#include <list>
#include <unordered_map>

#include <htslib/hfile.h>

#include "error.h"
#include "json.h"
#include "url.h"
#include "atom.h"
#include "sequence.h"

using std::set;
using std::copy;
using std::hash;
using std::list;
using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::string;
using std::vector;
using std::ostream;
using std::ifstream;
using std::ios_base;
using std::exception;
using std::to_string;
using std::make_pair;
using std::setprecision;
using std::unordered_map;

const ssize_t PEEK_BUFFER_CAPACITY(4096);
const int DEFAULT_FEED_CAPACITY(60);
const int DEFAULT_FEED_RESOLUTION(60);

class FeedSpecification {
friend ostream& operator<<(ostream& o, const FeedSpecification& specification);

public:
    const IoDirection direction;
    const int32_t index;
    URL url;
    Platform platform;
    int capacity;
    int resolution;
    uint8_t phred_offset;
    unordered_map< string, const HeadPGAtom > program_by_id;
    unordered_map< string, const HeadRGAtom > read_group_by_id;
    hFILE* hfile;

    FeedSpecification (
        const IoDirection& direction,
        const int32_t& index,
        const URL& url,
        const Platform& platform,
        const uint8_t& phred_offset);
    void set_capacity(const int& capacity);
    void set_resolution(const int& resolution);
    void register_rg(const HeadRGAtom& rg);
    void register_pg(const HeadPGAtom& pg);
    void describe(ostream& o) const;
    void probe();
};
ostream& operator<<(ostream& o, const FeedSpecification& specification);

class InputSpecification {
friend ostream& operator<<(ostream& o, const InputSpecification& specification);

public:
    Decoder decoder;
    bool disable_quality_control;
    bool long_read;
    bool include_filtered;
    list< URL > url_by_segment;
    vector< FeedSpecification* > feed_specification_by_segment;

    InputSpecification();
    void encode(Document& document, Value& node) const;
};

class ChannelSpecification {
friend ostream& operator<<(ostream& o, const ChannelSpecification& channel);

public:
    int32_t index;
    int64_t TC;
    kstring_t FS;
    kstring_t CO;
    Decoder decoder;
    bool disable_quality_control;
    bool long_read;
    bool include_filtered;
    bool undetermined;
    double concentration;
    Barcode multiplex_barcode;
    HeadRGAtom rg;
    list< URL > url_by_segment;
    vector< FeedSpecification* > feed_specification_by_segment;

    ChannelSpecification();
    ~ChannelSpecification();
    inline bool empty() const {
        return url_by_segment.empty();
    };
    string alias() const;
    void describe(ostream& o) const;
    void encode(Document& document, Value& node) const;
};
ostream& operator<<(ostream& o, const ChannelSpecification& specification);
void transcode_channel_specification(const Value& from, Value& to, Document& document);

#endif /* PHENIQS_MODEL_H */
