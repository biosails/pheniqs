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

#ifndef PHENIQS_MODEL_H
#define PHENIQS_MODEL_H

#include <set>
#include <unordered_map>

#include <htslib/hfile.h>
#include <htslib/kstring.h>

#include "error.h"
#include "json.h"
#include "constant.h"
#include "url.h"
#include "atom.h"
#include "sequence.h"

using std::set;
using std::copy;
using std::hash;
using std::setw;
using std::endl;
using std::cerr;
using std::cout;
using std::fixed;
using std::size_t;
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

class FeedSpecification {
friend ostream& operator<<(ostream& o, const FeedSpecification& specification);

public:
    const IoDirection direction;
    const size_t index;
    URL url;
    Platform platform;
    size_t capacity;
    size_t resolution;
    uint8_t phred_offset;
    unordered_map< string, const HeadPGAtom > program_by_id;
    unordered_map< string, const HeadRGAtom > read_group_by_id;
    hFILE* hfile;

    FeedSpecification (
        const IoDirection& direction,
        const size_t& index,
        const URL& url,
        const Platform& platform,
        const uint8_t& phred_offset);
    void set_capacity(const size_t& capacity);
    void set_resolution(const size_t& resolution);
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
    vector< URL > input_urls;
    vector< FeedSpecification* > feed_specifications;

    InputSpecification();
    void encode(Document& document, Value& node) const;
};

class ChannelSpecification {
friend ostream& operator<<(ostream& o, const ChannelSpecification& channel);

public:
    size_t index;
    size_t TC;
    kstring_t FS;
    kstring_t CO;
    Decoder decoder;
    bool disable_quality_control;
    bool long_read;
    bool include_filtered;
    bool undetermined;
    double concentration;
    Barcode multiplex_barcode;
    vector< URL > output_urls;
    HeadRGAtom rg;
    vector< FeedSpecification* > feed_specification;

    ChannelSpecification(size_t index);
    ~ChannelSpecification();
    inline bool writable() const {
        return output_urls.size() > 0;
    };
    string alias() const;
    void describe(ostream& o) const;
    void encode(Document& document, Value& node) const;
};
ostream& operator<<(ostream& o, const ChannelSpecification& specification);

#endif /* PHENIQS_MODEL_H */